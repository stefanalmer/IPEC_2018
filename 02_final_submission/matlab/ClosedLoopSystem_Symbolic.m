%%% symbolic expression of closed-loop system. Determine 
% coefficients of denominator

clear
syms s a b w1 w3 L1 L3
G = b/(s^2 + a*s + b);
H1 = L1*s/(s^2 + w1^2);
H3 = L3*s/(s^2 + w3^2);
H  = H1 + H3;
T = G*(1+H)/(1+G*H);

[n, d] = numden(T);
C = coeffs(d, s)