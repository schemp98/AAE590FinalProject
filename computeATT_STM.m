function A = computeATT_STM(q,w_c,I_T,I_C,h_wc)

Dq = D(q);
% Eqn (50)
CO = -skew(Dq*inv(I_c)*inv(Dq)*w_t);
FO = inv(I_T)*skew(h_wc)*I_T*Dq*inv(I_C)*inv(Dq);


A = zeros(6);
A(1:3,1:3) = -skew(w_c);
A(1:3,4:6) = q(4)*inv(Dq);
A(4:6,4:6) = FO + CO;