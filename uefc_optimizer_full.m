% Declare global variables here
% All constants should be globals
global m_pay rho g E min_SM C_mw max_elev_deflection rho_caps


m_pay = 0.3;
rho = 1.225;
g = 9.8066;
min_SM = 0.05;
C_mw = -0.13;
max_elev_deflection = 10*(pi/180);
rho_caps = 1000; % Find actual value for this


% Signatures can be worked out later

% The design is encoded as a (12,) array of design parameters:
% 1 - b_w,      Main wing span
% 2 - c_w,      Main wing mean chord
% 3 - Cl_nom,   Main wing nominal lift coefficient
% 4 - Cl_trim,  Main wing trim lift coefficient
% 5 - C_tw,     Main wing cap thickness
% 6 - C_ww,     Main wing cap width
% 7 - N,        Nominal load factor
% 8 - b_h,      Horizontal stab span
% 9 - c_h,      Horizontal stab mean chord
% 10 - Cl_hnom, Horizontal stab nominal lift coefficient
% 11 - x_h,     Horizontal stab NP abscissa from main wing NP abscissa
% 12 - SM_trim  Trim static margin

%Calculate CG shift
function[delta_x_pay] = get_delta_x_pay(x)

global min_SM C_mw
b_w = x(1);
c_w = x(2);
Cl_nom = x(3);
Cl_trim = x(4);
C_tw = x(5);
C_ww = x(6);
N = x(7);
b_h = x(8);
c_h = x(9);
Cl_hnom = x(10);
x_h = x(11);
SM_trim = x(12);

S_w = b_w*c_w;
S_h = b_h*c_h;
a_w = 2*pi/(1+(2/(b_w/c_w)));
a_h = 2*pi/(1+(2/(b_h/c_h)));
M_w = S_w*c_w*C_mw;

x_np = (x_h+S_h*a_h)/(S_h*a_h+S_w*a_w);
x_cg = (x_h*S_h*Cl_nom-M_w)/(Cl_nom*S_w+S_hC_lhnom);

delta_x_pay = x_np - c_w*SM_trim-x_cg;


end

%Calculate elevator trim constraint
function[con_elev_deflection] = get_elev_deflection(x)

global min_SM C_mw max_elev_deflection
b_w = x(1);
c_w = x(2);
Cl_nom = x(3);
Cl_trim = x(4);
C_tw = x(5);
C_ww = x(6);
N = x(7);
b_h = x(8);
c_h = x(9);
Cl_hnom = x(10);
x_h = x(11);
SM_trim = x(12);

S_w = b_w*c_w;
S_h = b_h*c_h;
a_w = 2*pi/(1+(2/(b_w/c_w)));
a_h = 2*pi/(1+(2/(b_h/c_h)));
M_w = S_w*c_w*C_mw;

x_np = (x_h+S_h*a_h)/(S_h*a_h+S_w*a_w);
x_cgtrim = x_np-c*SM_trim;
Cl_htrim = (-M_w-S_w*Cl_trim*x_cgtrim)/(S*h*(x_cgtrim-x_h));
delta_alpha = (Cl_trim-Cl_nom)/a_w;

elev_deflection = (Cl_htrim-Cl_hnom-delta_alpha*a*h)/a_h;
con_elev_deflection = elev_deflection - max_elev_deflection;

end

% Calculates velocity
function [v] = get_v(x)
global m_pay rho g
    b_w = x(1);
    c_w = x(2);
    N =  x(7);
    Cl_nom = x(3);

    W = get_m_tot(x) * g;
    S = b_w*c_w;
    v = sqrt((N*W)/((1/2)*rho*Cl_nom));
end

% Calculates drag force
function [F_d] = get_F_d(x,v)

end

% Calculate total aircraft mass
function [m_tot] = get_m_tot(x)
global g m_pay

    W_wing_tail_weight = get_wh(x);
    W_fusl = get_wfusl(x);
    W_pay = m_pay*g;

    m_tot = (W_wing_tail_weight + W_pay + W_fusl)/g;

end

% Calculate wing and tail weight
function [W_wing_tail_weight] = get_wh(x)
global g rho_foam
    Afac = 0.66;
    b_w = x(1);
    c_w = x(2);
    lam_w = x(3);
    tau_w = x(5);
    AR_w = b_w/c_w;
    S_w = b_w*c_w;
    b_h = x(8);
    c_h = x(9);
    lam_h = x(10);
    tau_h = x(12);
    AR_h = b_h/c_h;
    S_h = b_h*c_h;


    Wwing = (4/3)*Afac*rho_foam*g*tau_w*S_w^1.5*AR_w^(-0.5)*(lam_w^2+lam_w+1)/(lam_w+1)^2;
    Wtail = (4/3)*Afac*rho_foam*g*tau_h*S_h^1.5*AR_h^(-0.5)*(lam_h^2+lam_h+1)/(lam_h+1)^2;
    W_caps = 2*C_tw*C_ww*b_w*rho_caps*g;
    W_wing_tail_weight = Wwing + Wtail+W_caps;

end

% Calculate the weight of the fuselage
function [W_fusl] = get_wfusl(x)
global g
    mfuse0 = .185;
    mfusel = .060;
    mfuseS = .026;
    SPV = 0.225;
    bPV = 1.5;

    b_w = x(1);
    c_w = x(2);
    S_w = b_w*c_w;
    b_h = x(8);
    c_h = x(9);
    S_h = b_w*c_w;


    W_fusl = (mfuse0 + mfusel * ((b_w+b_h)/bPV) + mfuseS * ((S_w+S_h)/SPV))*g;

end


% This calculates the maximum prop thrust as a function of velocity
function [T_max] = get_T_max(x)
global m_pay rho g
    ct0 = 0.2093;
    ct1 = -0.2484;
    ct2 = -0.1386;
    Tmax_static = 2;
    Rprop       = 0.1016;
    Aprop       = np.pi*Rprop^2;
    Omega  = np.sqrt(Tmax_static/(0.5*rho*Rprop^2*Aprop*ct0));
    Lambda = V/(Omega*Rprop);
    CT = ct0+ct1*Lambda+ct2*Lambda^2;
    T_max = CT*0.5*rho*((Omega*Rprop)^2)*Aprop;
end

% Calculate the delta in cg position
function [delta_cg] = get_delta_cg(x)

end

% Calculate tip deflection
function [d_b] = get_d_b(x,m_tot)
global m_pay rho g
    b_w = x(1);
    c_w = x(2);
    tau_w = x(3);
    lam_w = x(4);
    tau_w = x(5);
    Cl_w = x(6);
    Ct_w = x(7);
    Cw_w = x(8);
    b_h = x(9);
    c_h = x(10);
    tau_h = x(11);
    lam_h = x(12);
    Cl_h = x(13);
    Ct_h = x(14);
    Cw_h = x(15);
    N = x(16);
    
    v = get_v(x);
    A = b_w * c_w;
    L_0 = (1/2)*Cl_w *(v^2)* rho * A;

    E = 3000000000;
    Gamma =  @(y) 1+((lam_w-1)/(lam_w + 1))-(((2.*y)/(b/2)).*((lam_w-1)/(lam_w+1)));
    I =@(y) (1/2).*Gamma(y).*((tau_w.*c_w)^2).*Ct_w.*Cw_w;
    Mx = @(y) (L_0/24).*(2.*b.^2.*sqrt(1-(4/b.^2).*y.^2) + 4.*(y.^2).*sqrt(1-(4/b.^2).*y.^2) - 3.*pi.*b.*y + 6.*b.*y.*asin(2.*y/b));
    u_doubleprime = @(y) arrayfun(@(y) (Mx(y))/(E.*I(y)),y);
    u_prime = @(y) arrayfun(@(y)integral(u_doubleprime, 0, y),y);
    u = integral(u_prime, 0, b/2);
    % y = b/2; Is this intentional? unused currently
    d_b = (eval(u)/b);

end

% Calculate turn radius
function [r_turn] = get_r_turn(x,v)
global m_pay rho g
    b_w = x(1);
    c_w = x(2);
    tau_w = x(3);
    lam_w = x(4);
    tau_w = x(5);
    Cl_w = x(6);
    Ct_w = x(7);
    Cw_w = x(8);
    b_h = x(9);
    c_h = x(10);
    tau_h = x(11);
    lam_h = x(12);
    Cl_h = x(13);
    Ct_h = x(14);
    Cw_h = x(15);
    N = x(16);

    v = get_v(x)
    r_turn = ((v^2)*sqrt((N^2) - 1))/g;;

end

% Calculate the complete objective
% Remember to normalize, we can get approximate scaling factors by
% optimizing for each objective separately
function [obj] = get_obj(x)

end

% We will likely use the ga optimizer.
% Due to the way it works, we have to report the optimizer constaints as a
% 1D vector in ceq.
% c should not be used unless we need equality constraints
function [c,ceq] = get_constraints(x)

c = [];

end

% Optimization bounds
print = get_d_b([1,1,1,1,1,1,1,1,1,1,1,1])
% Optimizer run parameters
