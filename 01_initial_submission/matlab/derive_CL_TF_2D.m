

% ---------------------------------
% --- 2D problem ------------------
% ---------------------------------

s = sym('s');
a = sym('a');
b = sym('b');
w1 = sym('w1');
w3 = sym('w3');
L1 = sym('L1');
L3 = sym('L3');

G = b/(s^2 + a*s + b);
H = L1*s/(s^2 + w1^2) + L3*s/(s^2 + w3^2);
Gcl = simple( G*(1+H)/(1+G*H) );
% redult: (b*s^4 + (L1*b + L3*b)*s^3 + (b*w1^2 + b*w3^2)*s^2 + (L3*b*w1^2 + L1*b*w3^2)*s + b*w1^2*w3^2)/(s^6 + a*s^5 + (w1^2 + w3^2 + b)*s^4 + (a*w1^2 + a*w3^2 + L1*b + L3*b)*s^3 + (w1^2*w3^2 + b*w1^2 + b*w3^2)*s^2 + (a*w1^2*w3^2 + L3*b*w1^2 + L1*b*w3^2)*s + b*w1^2*w3^2)
 
% explicit expression for denominator
% den = [1 ;
% a ;
% (w1^2 + w3^2 + b);
% (a*w1^2 + a*w3^2 + L1*b + L3*b);
% (w1^2*w3^2 + b*w1^2 + b*w3^2);
% (a*w1^2*w3^2 + L3*b*w1^2 + L1*b*w3^2) ;
% b*w1^2*w3^2];
