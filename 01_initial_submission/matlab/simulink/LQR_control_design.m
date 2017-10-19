       

% local variable names
Tc = VSI.Treg;  % sampling time of controller
wg = 2*pi*fn;   % grid frequency
L = 85e-6;      % filter inductance 
r = VSI.Rf;     % inductor parasitic resistance
C = 275e-6;     % filter capacitance
rc = VSI.Rcf;   % capacitor parasitic resistance

CTRL.r = r;
CTRL.wg = wg;
CTRL.L = L;
CTRL.C = C;


%--- system dynamics -----------------------------------

% let system state be in abc-coordinates. LC filter state is defiend as
% x = [ia ib ic va vb vc] where ix are the filter currents and vx are the
% filter voltage.
% Let w = [wa wb wc] be load current and vDC be dc-link voltage
% and u = [ua ub uc] be inverter voltage
% system dynamics are
% \dot x = Ac*x + Bc*vDC*s(t) + Fc*w
% where system matrices are defined below
%
% Note: consider case with common connection point of filter capacitors 
% is connected to neutral point of DC link. 
Ac = [-1/L*(r+rc), -1/L;
       1/C,  0];
Bc = [1/L; 0];
Fc = [rc/L; -1/C];
% -----------------------------------------------------------


% --- harmonic balance --------------------------------------
% define matrices (Au,Bu) and (Ai,Bi) such that
% u_ref = Au*v_ref Bu*w
% i_ref = Ai*v_ref Bi*w
% where v_ref = V*[cos(wt), sin(wt)] is the (given)
% reference for the capacitor voltage and
% w= W*[cos(wt), sin(wt)] is the (given) load current 
% and i_ref = I*[cos(wt), sin(wt)] is the inductor current reference
% and u_ref = U*[cos(wt), sin(wt)] is hte inverter voltage 
%
Jss = [0 -1; 1 0];
I = eye(2);

CTRL.Au = -wg^2*L*C*I + r*wg*C*Jss + I;
CTRL.Bu = wg*L*Jss + r*I;

CTRL.Ai = wg*C*Jss;
CTRL.Bi = I;
% -----------------------------------------------------------


% --- LQR feedback design ----------------------------------
% discrete time model
Ad = expm(Ac*Tc);
Bd = [eye(2),zeros(2,1)]*expm([Ac,Bc;zeros(1,3)]*Tc)*[zeros(2,1);1];

Qlq = zeros(2);
Qlq(1,1) = 1; 
Qlq(2,2) = 4; 
Rlq = 17; 
[Klq,Slq,Elq] = dlqr(Ad,Bd,Qlq,Rlq);
CTRL.Klq = Klq;
% -----------------------------------------------------------


% --- proportional resonant controller for offset free tracking
% continuous transfer functions:
% lam*s/(s^2 + (n*w)^2)

%lam = 611.8312;
% optimal gains from 4D-optimization 
% 111.1867
% 313.8555
% 344.3558
% 690.0204

  
%lam = 20;
%lam1_pr = 150.4338;
%lam3_pr = 428.4452;
%lam5_pr = 0;
%lam7_pr = 0;
%lam9_pr = 0;

%lam1_pr = lam;
%lam3_pr = lam;
%lam5_pr = lam;
%lam7_pr = lam;
%lam9_pr = lam;

w1_pr = 50*2*pi;
w3_pr = 3*50*2*pi;
w5_pr = 5*50*2*pi;
w7_pr = 7*50*2*pi;
w9_pr = 9*50*2*pi;

s=tf('s');
PR_1 = lam1_pr*s/(s^2 + w1_pr^2);
PR_3 = lam3_pr*s/(s^2 + w3_pr^2);
PR_5 = lam5_pr*s/(s^2 + w5_pr^2);
PR_7 = lam7_pr*s/(s^2 + w7_pr^2);
PR_9 = lam9_pr*s/(s^2 + w9_pr^2);

% discretize time representation
PR_1d = c2d(PR_1,Tc,'zoh');
PR_3d = c2d(PR_3,Tc,'zoh');
PR_5d = c2d(PR_5,Tc,'zoh');
PR_7d = c2d(PR_7,Tc,'zoh');
PR_9d = c2d(PR_9,Tc,'zoh');

[num_1 den_1] = tfdata(PR_1d);
CTRL.num_1 = num_1{1};
CTRL.den_1 = den_1{1};
[num_3 den_3] = tfdata(PR_3d);
CTRL.num_3 = num_3{1};
CTRL.den_3 = den_3{1};
[num_5 den_5] = tfdata(PR_5d);
CTRL.num_5 = num_5{1};
CTRL.den_5 = den_5{1};
[num_7 den_7] = tfdata(PR_7d);
CTRL.num_7 = num_7{1};
CTRL.den_7 = den_7{1};
[num_9 den_9] = tfdata(PR_9d);
CTRL.num_9 = num_9{1};
CTRL.den_9 = den_9{1};
% -----------------------------------------------------------


% --- low pass filter for the load current estimation -------
w_lp = 1000*2*pi; % 1 kHz low pass filter
s=tf('s');
P_lp = w_lp/(s + w_lp);

% discretize time representation
P_lpd = c2d(P_lp,Tc,'zoh');
[num_lp den_lp] = tfdata(P_lpd);
CTRL.num_lp = num_lp{1};
CTRL.den_lp = den_lp{1};
% -----------------------------------------------------------


% --- Kalman fitler estimator for load current --------------
% model for the load current harmonics
% need following model structure for Kalman filter design
% x = Ax + Bu + Gww            {State equation}
% y = Cx + Du + Hww + vv        {Measurements}

% KEST uses [u(t);y(t)] to generate optimal estimates y_e(t),x_e(t) of 
% y(t),x(t) by:
%  .
%  x_e  = Ax_e + Bu + L (y - Cx_e - Du)
% 
%  |y_e| = | C | x_e + | D | u
%  |x_e|   | I |       | 0 | 

A1 = [0,-1*wg; 1*wg,0];
BB = [0;0];
GG = eye(2);
CC = [1, 0];
DD = 0;
HH = [0,0];
sys1 = ss(A1,[BB GG],CC,[DD HH]);

Qk = 1e4*eye(2);
Rk = 1e-6*eye(1); 
NN = [0; 0];

% kalman filter synthesis
[KEST1, L1, P1] = kalman(sys1,Qk,Rk,NN);

% estimator realizations
AK1 = A1 - L1*CC;
BK1 = L1;
CK1 = eye(2);
DK1 = [0; 0];

% discrete time Kalman fitler
% discretize the continuous time estimator
CTRL.AKd1 = expm( (A1 - L1*CC)*Tc );
CTRL.BKd1 = [eye(2),zeros(2,1)]*expm( [(A1 - L1*CC),L1;zeros(1,3)]*Tc)*[zeros(2,1);1];
CTRL.CKd1 = eye(2);
CTRL.DKd1 = [0; 0];
% ------------------------------------------------------------------

% Added 20170626 ---------------------------------------------------
% Kalman filter design in discrete  time

A1 = [0,-1*wg; 1*wg,0];
BB = [0;0];
GG = eye(2);
CC = [1, 0];
DD = 0;
HH = [0,0];
%sys1 = ss(A1,[BB GG],CC,[DD HH]);
sys1_d = c2d(sys1,Tc);

Qk = 1e4*eye(2);
Rk = 1e-6*eye(1); 
NN = [0; 0];

[KEST1_d, L1_d, P1_d] = kalman(sys1_d,Qk,Rk,NN,'current');
CTRL.AKd2 = KEST1_d.a;
CTRL.BKd2 = KEST1_d.b(:,2);
CTRL.CKd2 = KEST1_d.c(2:3,:);
CTRL.DKd2 = KEST1_d.d(2:3,2);

% -----------------------------------------------------------------
