%%% generate the Bode plot of the closed-loop system with 2 PR ctrller

% Closed-loop system - pole location and bode plot
L    =  85e-6;   % [H]
Ccap = 275e-6;   % [F]
R    = 0.010;    % [ohm]

% G = b/(s^2 + a*s + b)
b = 10463/(Ccap*10000*L);
a = (3259 + 10000*R)/(10000*L);

s = tf('s');
G = b/(s^2 + a*s + b);

% PR controller
w1 = 2*pi*50;
w3 = 3*w1;

L1 = 500;
L3 = 200;

H1 = L1*s/(s^2 - w1^2);
H3 = L3*s/(s^2 - w3^2);

T = G*(1-H1-H3)/(1-G*(H1+H3));

pole(T)
zero(T)
pzmap(T)

figure
bode(T)
