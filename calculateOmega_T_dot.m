function wDot = calculateOmega_T_dot(q,w_t,I_c,T_cc,h_wc)

Dq = D(q);
% Equation (50)
CO = -skew(Dq*inv(I_c)*inv(Dq)*w_t);

% Equation (21)
term1 = CO*I_c*inv(Dq)*w_t;
term2 = CO*h_wc;
term3 = Dq*inv(I_c)*T_cc;

wDot  = term1 + term2 + term3;

end