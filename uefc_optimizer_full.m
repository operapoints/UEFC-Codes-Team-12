warning('off', 'all');
% Declare global variables here
% All constants should be globals
global m_pay rho g min_SM C_mw max_elev_deflection rho_caps tau lam rho_balsa mu rho_foam t_h spaneff mass_margin drag_margin


m_pay = 0.3;
rho = 1.225;
g = 9.8066;
min_SM = 0.05;
C_mw = -0.13;
max_elev_deflection = 10*(pi/180);
rho_caps = 500; 
tau = 0.12;
lam = 0.5;
rho_balsa = 500;
rho_foam = 32;
mu = 1.8e-5;
t_h = 1.6e-3;
spaneff = 0.95;
mass_margin = 1.1;
drag_margin = 1.1;


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
% 13 - N_trim   Trim load factor

%Calculate CG shift
function[delta_x_pay] = get_delta_x_pay(x, m_tot)

global min_SM C_mw g m_pay
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

W_tot = m_tot * g;

S_w = b_w*c_w;
S_h = b_h*c_h;
a_w = 2*pi/(1+(2/(b_w/c_w)));
a_h = 2*pi/(1+(2/(b_h/c_h)));
M_w = S_w*c_w*C_mw;

x_np = (x_h*S_h*a_h)/(S_h*a_h+S_w*a_w);
x_cg = (x_h*S_h*Cl_hnom-M_w)/(Cl_nom*S_w+S_h*Cl_hnom);

delta_x_pay = (x_np - (c_w*SM_trim+x_cg))*(m_tot/m_pay);


end

%Calculate elevator trim constraint
function[con_elev_deflection,Cl_htrim] = get_elev_deflection(x)

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

x_np = (x_h*S_h*a_h)/(S_h*a_h+S_w*a_w);
x_cgtrim = x_np-c_w*SM_trim;
Cl_htrim = (-M_w-S_w*Cl_trim*x_cgtrim)/(S_h*(x_cgtrim-x_h));
delta_alpha = (Cl_trim-Cl_nom)/a_w;

elev_deflection = (Cl_htrim-Cl_hnom-delta_alpha*a_h)/a_h;
con_elev_deflection = elev_deflection - max_elev_deflection;

end

% Calculates velocity
function [v] = get_v(x,m_tot)
global m_pay rho g
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

    W = m_tot * g;
    S_w = b_w*c_w;
    S_h = b_h*c_h;
    v = sqrt((2*N*W)/(rho*(S_w*Cl_nom+S_h*Cl_hnom)));
end

% Calculates drag force
function [F_d] = get_F_d(x,v)
    global rho mu tau spaneff mass_margin
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
    q = 0.5*rho*v^2;
    S_w = c_w*b_w;
    S_h = c_h*b_h;

    SPV       = 0.225 ;
    CDA_fuse0 = 0.002;
    CDA_fuseS = 0.002;
    CDfuse = (1/S_w) * (CDA_fuse0 + CDA_fuseS*(S_w/SPV));
    F_d_fuse = q*S_w*CDfuse;
    
    Re_w = rho*v*c_w/mu;
    cd0    = 0.020*(1+tau^2);
    cd1    = -0.005/(1+6*tau);
    cd2    = 0.160/(1+60*tau);
    cd8    = 1.0;
    cl0    = 1.25 - 3*tau;
    Re_ref = 1E5;
    Re_a   = -0.75;
    cl2d   = Cl_nom;
    cdpfac = cd0 + cd1*(cl2d-cl0) + cd2*(cl2d-cl0)^2 + cd8*(cl2d-cl0)^8;
    CDp_w    = (cdpfac*(Re_w/Re_ref)^Re_a);
    F_dp_w = q*S_w*CDp_w;

    CDi_w = (Cl_nom^2)/(pi*spaneff*(b_w/c_w));
    F_di_w = q*S_w*CDi_w;


    CDp_h    = 0.025;
    F_dp_h = 1.5*q*S_h*CDp_h;

    CDi_h = (Cl_hnom^2)/(pi*spaneff*(b_h/c_h));
    F_di_h = q*S_h*CDi_h; % Factor of 1.5 to account for the rudder

    F_d = (F_d_fuse+F_di_w+F_dp_w+F_di_h+F_dp_h)*mass_margin;

end

% Calculates trim drag force
function [F_dtrim] = get_F_dtrim(x,v,Cl_htrim)
    global rho mu tau spaneff mass_margin
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
    q = 0.5*rho*v^2;
    S_w = c_w*b_w;
    S_h = c_h*b_h;

    SPV       = 0.225 ;
    CDA_fuse0 = 0.002;
    CDA_fuseS = 0.002;
    CDfuse = (1/S_w) * (CDA_fuse0 + CDA_fuseS*(S_w/SPV));
    F_d_fuse = q*S_w*CDfuse;
    
    Re_w = rho*v*c_w/mu;
    cd0    = 0.020*(1+tau^2);
    cd1    = -0.005/(1+6*tau);
    cd2    = 0.160/(1+60*tau);
    cd8    = 1.0;
    cl0    = 1.25 - 3*tau;
    Re_ref = 1E5;
    Re_a   = -0.75;
    cl2d   = Cl_nom;
    cdpfac = cd0 + cd1*(cl2d-cl0) + cd2*(cl2d-cl0)^2 + cd8*(cl2d-cl0)^8;
    CDp_w    = (cdpfac*(Re_w/Re_ref)^Re_a);
    F_dp_w = q*S_w*CDp_w;

    CDi_w = (Cl_trim^2)/(pi*spaneff*(b_w/c_w));
    F_di_w = q*S_w*CDi_w;


    CDp_h    = 0.025;
    F_dp_h = 1.5*q*S_h*CDp_h;

    CDi_h = (Cl_htrim^2)/(pi*spaneff*(b_h/c_h));
    F_di_h = q*S_h*CDi_h; % Factor of 1.5 to account for the rudder

    F_dtrim = (F_d_fuse+F_di_w+F_dp_w+F_di_h+F_dp_h)*mass_margin;

end

% Calculate total aircraft mass
function [m_tot] = get_m_tot(x)
global g m_pay mass_margin

    W_wing_tail_weight = get_wh(x);
    W_fusl = get_wfusl(x);
    W_pay = m_pay*g;

    m_tot = (mass_margin*(W_wing_tail_weight + W_pay + W_fusl))/g;

end

% Calculate wing and tail weight
function [W_wing_tail_weight] = get_wh(x)
global g rho_foam lam tau rho_balsa t_h

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

    Afac = 0.66;

    lam_w = lam;
    tau_w = tau;
    AR_w = b_w/c_w;
    S_w = b_w*c_w;

    S_h = b_h*c_h;


    Wwing = (4/3)*Afac*rho_foam*g*tau_w*S_w^1.5*AR_w^(-0.5)*(lam_w^2+lam_w+1)/(lam_w+1)^2;
    W_caps = 2*C_tw*C_ww*b_w*rho_balsa*g;
    W_tail = (4/3)*Afac*rho_foam*g*tau_w*S_h^1.5*(b_h/c_h)^(-0.5)*(lam_w^2+lam_w+1)/(lam_w+1)^2;
    W_caps_h = 2*C_tw*C_ww*b_h*rho_balsa*g;
    W_wing_tail_weight = Wwing + W_caps + 1.5*W_tail+1.5*W_caps_h;

end

% Calculate the weight of the fuselage
function [W_fusl] = get_wfusl(x)
global g
    mfuse0 = .245;
    mfusel = .060;
    mfuseS = .026;
    SPV = 0.225;
    bPV = 1.5;

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

    fuse_len_delta = (x_h - 1)*2;
    fuse_len_delta_m = 0.0425*fuse_len_delta;
    


    W_fusl = (mfuse0 + mfusel * ((b_w)/bPV) + mfuseS * ((S_w)/SPV)+fuse_len_delta_m)*g;

end


% This calculates the maximum prop thrust as a function of velocity
function [T_max] = get_T_max(v)
global m_pay rho g
    ct0 = 0.2093;
    ct1 = -0.2484;
    ct2 = -0.1386;
    Tmax_static = 2;
    Rprop       = 0.1016;
    Aprop       = pi*Rprop^2;
    Omega  = sqrt(Tmax_static/(0.5*rho*Rprop^2*Aprop*ct0));
    Lambda = v/(Omega*Rprop);
    CT = ct0+ct1*Lambda+ct2*Lambda^2;
    T_max = CT*0.5*rho*((Omega*Rprop)^2)*Aprop;
end

% Calculate tip deflection
function [d_b] = get_d_b(x,m_tot)
global m_pay rho g lam tau
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
    

    L_0 = m_tot*((b_w*c_w*Cl_nom)/(Cl_nom*b_w*c_w+Cl_hnom*b_h*c_h)) * g * N * (4/(pi*b_w));

    E = 3e+9;
    Gamma =  @(y) 1+((lam-1)/(lam + 1))-(((2.*y)/(b_w/2)).*((lam-1)/(lam+1)));
    I =@(y) 2*(((1/2).*Gamma(y).*tau.*c_w-(C_tw/2))^2).*C_tw.*C_ww;
    Mx = @(y) (L_0/24).*(2.*b_w.^2.*sqrt(1-(4/b_w.^2).*y.^2) + 4.*(y.^2).*sqrt(1-(4/b_w.^2).*y.^2) - 3.*pi.*b_w.*y + 6.*b_w.*y.*asin(2.*y/b_w));
    u_doubleprime = @(y) arrayfun(@(y) (Mx(y))/(E.*I(y)),y);
    u_prime = @(y) arrayfun(@(y)integral(u_doubleprime, 0, y),y);
    u = integral(u_prime, 0, b_w/2);
    % y = b/2; Is this intentional? unused currently
    d_b = (u/b_w);

end

% Calculate tip deflection for tail
function [d_b_h] = get_d_b_h(x,m_tot,Cl_htrim)
global m_pay rho g lam tau
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
    N_trim = x(13);
    
    %The factor of 1.5 in front is to account for control forces
    L_0a = 1*m_tot*((b_h*c_h*Cl_htrim)/(Cl_trim*b_w*c_w+Cl_htrim*b_h*c_h)) * g * N_trim * (4/(pi*b_h));
    L_0b = 1*m_tot*((b_h*c_h*Cl_hnom)/(Cl_nom*b_w*c_w+Cl_hnom*b_h*c_h)) * g * N * (4/(pi*b_h));
    L_0 = max(abs(L_0a),abs(L_0b));
    %_h = c_h
    %meanhthickness = 0.5*tau*c_h-(C_tw/2)
    E = 3e+9;
    Gamma =  @(y) 1+((lam-1)/(lam + 1))-(((2.*y)/(b_h/2)).*((lam-1)/(lam+1)));
    I =@(y) 2*((0.5.*Gamma(y).*tau.*c_h-(C_tw/2))^2).*C_tw.*C_ww;
    Mx = @(y) (L_0/24).*(2.*b_h.^2.*sqrt(1-(4/b_h.^2).*y.^2) + 4.*(y.^2).*sqrt(1-(4/b_h.^2).*y.^2) - 3.*pi.*b_h.*y + 6.*b_h.*y.*asin(2.*y/b_h));
    u_doubleprime = @(y) arrayfun(@(y) (Mx(y))/(E.*I(y)),y);
    u_prime = @(y) arrayfun(@(y)integral(u_doubleprime, 0, y),y);
    u = integral(u_prime, 0, b_h/2);
    % y = b/2; Is this intentional? unused currently
    d_b_h = (u/b_h);

end

% Calculate turn radius
function [r_turn] = get_r_turn(x,v)
global m_pay rho g
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

    r_turn = ((v^2)/(g*sqrt((N^2) - 1)));

end

% Calculate the complete objective
% Remember to normalize, we can get approximate scaling factors by
% optimizing for each objective separately
function [obj] = get_obj(x)
try

    m_tot = get_m_tot(x);
    v = get_v(x,m_tot);
    delta_x_pay = get_delta_x_pay(x,m_tot);
    
    %obj = -v;
    %obj = -delta_x_pay;
    obj = -(v/8.9987+delta_x_pay/0.7918);
    % TODO: Normalize obj by sub objectives
catch
    obj = 1e6;
end
end

% We will likely use the ga optimizer.
% Due to the way it works, we have to report the optimizer constaints as a
% 1D vector in ceq.
% c should not be used unless we need equality constraints
function [c,ceq] = get_constraints(x)
try    
    global g rho tau
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
    N_trim = x(13);

    m_tot = get_m_tot(x);
    v = get_v(x,m_tot);
    T_max = get_T_max(v);
    F_d = get_F_d(x,v);
    d_b = get_d_b(x,m_tot);
    con_thrust_drag = F_d - T_max;
    con_d_b = d_b-0.05;
    r_turn = get_r_turn(x,v);
    con_r_turn = r_turn - 12.5;
    [con_elev_deflection, Cl_htrim] = get_elev_deflection(x);
    W = m_tot*g;
    v_trim = sqrt((2*N_trim*W)/(rho*(b_w*c_w*Cl_trim+b_h*c_h*Cl_htrim)));
    T_trim = get_T_max(v_trim);
    F_dtrim = get_F_dtrim(x,v_trim,Cl_htrim);
    con_trim_thrust = F_dtrim - T_trim;
    r_turn_trim = ((v_trim^2)/(g*sqrt((N_trim^2) - 1)));
    con_trim_r_turn = r_turn_trim - 12.5;
    d_b_h = get_d_b_h(x,m_tot,Cl_htrim);
    con_d_b_h = d_b_h-0.02;
    con_h_thickness = 2*C_tw - tau*c_h;



    
   intm = [con_thrust_drag;
        con_d_b;
        con_r_turn;
        con_elev_deflection;%Zeroed for delta_xpay convergence
        con_trim_thrust;%Zeroed for delta_xpay convergence
        con_trim_r_turn;
        con_d_b_h;
        con_h_thickness;
        ];


    has_invalid = any(isnan(intm) | isinf(intm));
    if has_invalid
        
        ceq = [];
        c = [1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6];
    else
        ceq = [];
        c = intm;
    end
catch
    ceq = [];
    c = [1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6];
end

end



%    b_w = x(1);
%    c_w = x(2);
%    Cl_nom = x(3);
%    Cl_trim = x(4);
%    C_tw = x(5);
%    C_ww = x(6);
%    N = x(7);
%    b_h = x(8);
%    c_h = x(9);
%    Cl_hnom = x(10);
%    x_h = x(11);
%    SM_trim = x(12);
%    N_trim = x(13);

% Optimization bounds
%[print1,print2] = get_constraints_debug([1.1,0.13,0.7,0.65,0.002,0.005,1.8,0.15,0.05,0,1,0.05])
options = optimoptions('fmincon', 'Display', 'iter','Algorithm','interior-point','ConstraintTolerance',1e-6);
%options = optimoptions('ga', 'Display', 'iter');
%options.TolCon = 0.03;
%options.PopulationSize = 20;
intcon = [];


lb = [0,0,0,0,0,0,1,0,0,-0.6,1,0.05,1];
ub = [4,0.5,0.8,0.8,0.003,0.01,1.5,5,0.5,0.8,1.5,0.5,1.5];
%[x,opt]=ga(@get_obj,13,[],[],[],[],lb,ub,@get_constraints,intcon,options)
x0 = [2.2435;
    0.2005;
    0.8000;
    0.8000;
    0.0017;
    0.0100;
    1.0935;
    1.3226;
    0.0967;
    0.2043;
    1.5000;
    0.0500;
    1.0659];
%diff = x0-ub;
[x,opt]=fmincon(@get_obj,x0,[],[],[],[],lb,ub,@get_constraints,options)

disp(x)
%x_dbg = x
%x_dbg = [2.0256,0.1236,0.8,0.6655,0.0050,0.0020,1.214,0.4,0.2,-0.1379,2,0.05,1.2019]
obj = get_obj(x)
[c,ceq]=get_constraints(x)


