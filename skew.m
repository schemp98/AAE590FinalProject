function S = skew(vec)
% https://en.wikipedia.org/wiki/Skew-symmetric_matrix#Cross_product

S       = zeros(3);
S(1,2)  = -vec(3); S(2,1)   = -S(1,2);
S(1,3)  =  vec(2); S(3,1)   = -S(1,3);
S(2,3)  = -vec(1); S(3,2)   = -S(2,3);