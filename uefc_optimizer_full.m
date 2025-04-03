% Declare global variables here
% All constants should be globals
global m_pay rho g

m_pay = 0.3;
rho = 1.225;
g = 9.8066;

% Signatures can be worked out later

% The design is encoded as a (15,) array of design parameters:
% b_w,c_w,tau_w,

% Calculates velocity
function [v] = get_v(x)

end

% Calculates drag force
function [F_d] = get_F_d(x,v)

end

% Calculate aircraft     mass
function [m_tot] = get_m_tot(x)

end

% This calculates the maximum prop thrust as a function of velocity
function [T_max] = get_T_max(x)

end

% Calculate the delta in cg position
function [delta_cg] = get_delta_cg(x)

end

% Calculate tip deflection
function [del_tip] = get_del_tip(x,m_tot)

end

% Calculate turn radius
function [r_turn] = get_r_turn(x,v)

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

% Optimizer run parameters

