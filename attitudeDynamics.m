function dx = attitudeDynamics(t,x0,I_c,T_cc,h_wc)

q0 = x0(1:4);
q0 = q0/norm(q0);  % Ensure that this stays normalized
w0 = x0(5:end);   % relative body rate in the Target Frame

qd = relQuatRate(t,q0,w0);  % Calculate Quaternion Rate

% T_cc = zeros(3,1);  %  SCC -Test until I can figure out how to calculate derivative of angular momentum of recation wheels (hd_wc)
wd = calculateOmega_T_dot(q0,w0,I_c,T_cc,h_wc);
dx = [qd;wd];