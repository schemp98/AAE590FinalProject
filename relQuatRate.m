function dq = relQuatRate(t,q,w_c)


% Eqn (47)
Omega =@(w_c) [ skew(w_c) w_c; -w_c' 0];


% Eqn (46)  - SCC - Not sure if this is needed...
dq = 0.5*Omega(w_c)*q;

