function dx = coupledDynamics(t,x0,mu, a, e, controller_type,B,R,Q,I_c,Pi_c,Pj_t)
%  State rho_ij
x0_rho_ij       = x0(1:3);
x0_rho_ij_d     = x0(4:6);
% x0_pos          = x0(1:6); % rho0; rho0d
x0_q            = x0(7:10);  % relative attitude
x0_w            = x0(11:13);  % rate
x0_att          = [x0_q;x0_w];

DCM             = D(x0_q);

use_RWA = true;
if  use_RWA

h_wc  = x0(14:16);
GO_factor = 1;
else
GO_factor = 1;
    h_wc  = x0(14:16)*0;
end


calc_Junction = (norm(Pi_c) + norm(Pj_t)) > 0;
calc_Att   =  norm([norm(x0_q(1:3))  norm(x0_w)]) > 1e-20;

if not(calc_Att)
shaun = 1;
end

if calc_Junction
Pi_c_t      = DCM*Pi_c;  % Point on Chaser in Target Frame;
rho0        = x0_rho_ij - Pi_c_t + Pj_t;   % eqn (26) 
wXP         = cross(x0_w,Pi_c_t);          
rhod0       = x0_rho_ij_d - wXP;            % eqn (27)
else
rho0 = x0(1:3);
rhod0 = x0(4:6);
end

x0_pos = [rho0;rhod0];


% calculate dynamics for relative Center of Masses
%  be sure to remove the control input from this function
dx_rho0 = linearized_rel_orbital_dynamics(t, x0_pos, mu, a, e); % rho0d; rho0dd
K = calculateControllerGainfunction(t,x0_pos, mu, a, e, controller_type,B,R,Q);
u_pos      = K*[x0(1:6)];

if calc_Att
% calculate dynamics for relative Attitude
I_t = I_c;  % Assuming Target and Chaser are the same shape
A_att = computeATT_STM(x0_q,x0_w,I_t,I_c,h_wc);

GO = GO_factor*DCM*inv(I_c);  % Think if we are using RWAs this is negative Eqn (52)
B_att = [zeros(3);GO];

Q_att = eye(6)*1000; R_att = 0.1*eye(3);

HM_att = [A_att -B_att*inv(R_att)*B_att';-Q_att -A_att'];
[Y,X] = ric_schr(HM_att);

P = X*inv(Y);
K = -inv(R_att)* B_att'* P;

% Use small angle approximation
phi   = DCM(2,3);  % roll
theta = DCM(1,3);  % pitch
psi   = DCM(2,1);  % yaw


x_att_K = [phi;theta;psi;x0_w];
T_cc = K*x_att_K;
else
T_cc = zeros(3,1);
end

h_wcd = -T_cc;

dx_att = attitudeDynamics(t,x0_att,I_c,T_cc,h_wc);
wd = dx_att(5:7);

% calculate dynamics for relative junction points
if calc_Junction
rho_ij_dd = dx_rho0(4:6) + cross(wd,Pi_c_t) + cross(x0_w,wXP) ; %eqn (28)
else
rho_ij_dd = dx_rho0(4:6);
end


dx = [x0_rho_ij_d;    
    rho_ij_dd + u_pos ;
    dx_att;
    h_wcd];