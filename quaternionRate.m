function qDot = quaternionRate(q,w)

% Eqn 12
Q =@(q) [eye(3)*q(4) + skew(q(1:3));q(1:3)'];


qDot = 0.5*Q(q)*w;