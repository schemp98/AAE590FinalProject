function A = computePOS_Far_STM(mu,r_c,w_OT,wD_OT)

Dq = D(q);
% Eqn (42)


A = zeros(6);
% A(1:3,1:3) = -skew(w_c);
A(1:3,4:6) = eye(3);

A(4:6,1:3) = skew(-[0;wD_OT;0]) + diag(ones(3,1)*-mu/r_c^3)...
    + diag([1;1;0]*w_OT^2);
A(4:6,4:6) = skew(-[0;w_OT;0]);
