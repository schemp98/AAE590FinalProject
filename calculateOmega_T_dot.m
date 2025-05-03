function wDot = calculateOmega_T_dot(q,w_t,I_c,T_cc,h_wc)

I_T = I_c;  % Assume that Target and Chaser are the same size
Dq = D(q);
% Equation (50)
CO = -skew(Dq*inv(I_c)*inv(Dq)*w_t)*I_c*inv(Dq);

% Equation (21)
term1 = CO*w_t;
% term2 = inv(I_T)*skew(h_wc)*I_T*Dq*inv(I_c)*inv(Dq);
term2 = -skew(Dq*inv(I_c)*inv(Dq)*w_t)*h_wc;
term3 = Dq*inv(I_c)*T_cc;

wDot  = term1 + term2 + term3;

end