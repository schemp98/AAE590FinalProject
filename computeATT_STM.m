function A = computeATT_STM(q,w_t,I_T,I_C,h_wc)

Dq = D(q);
w_c = Dq'*w_t; % Transform Relative Rate from Target Frame to Chaser Frame
% Eqn (50)
CO = -skew(Dq*inv(I_C)*inv(Dq)*w_t)*I_C*inv(Dq);
FO = inv(I_T)*skew(h_wc)*I_T*Dq*inv(I_C)*inv(Dq);


A = zeros(6);
A(1:3,1:3) = -skew(w_c);
A(1:3,4:6) = q(4)*inv(Dq);
A(4:6,4:6) = FO + CO;