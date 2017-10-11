%%% Derive a second-order model of a standard controlled converter
close all
clear 
clc

version = 2; % 1: PI, 2: P

%%% 1. cascaded PI controllers in abc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if version == 1
% Filter parameters
L    =  85e-6;   % [H]
Ccap = 275e-6;   % [F]
R    =  0.010;   % [ohm]

% PI controller parameters (kp + kp/Ti*s)
kpI = 0.427;
TiI = 0.000248;
kpU = 0.2765;
TiU = 0.0099;

%%% system equations
% dudt = 1/C*i
% didt = -1/L*u-R/L*i+1/L*u0
% dx1dt = iRef - i
% dx2dt = uRef - u
% iRef = kpU*(uRef-u)+kpU/TiU*x2
% u0 = kpI*(iRef-i) + kpI/TiI*x1 = -kpI*i + kpI/TiI*x1 + kpI*(kpU*(uRef-u)+kpU/TiU*x2)

%%% state: x = [u i x1 x2]
% dxdt = A*x+B*u
%    y = C*x

A = [ 0,            1/Ccap ,    0,          0;
    -1/L-kpI*kpU/L, -R/L-kpI/L, kpI/TiI/L,  kpI*kpU/TiU/L;
    -kpU,           -1,         0,          kpU/TiU;
    -1,             0,          0,          0];
    
 B = [ 0;           kpI*kpU/L ; kpU ;       1];
 
 C = [ 1,           0,          0,          0];
 
 sys = ss(A, B, C, 0);
 
 bode(sys), grid on
 
 pole(sys)
 
 figure
 pzmap(sys), axis equal
end
 
%%% 2. Cascaded P controllers in abc (PRs to be added) %%%%%%%%%%%%%%%%%%%
if version == 2
    disp(['INFO (' mfilename '.m): Considering cascaded P controllers'])
% Filter parameters
L    =  85e-6;   % [H]
Ccap = 275e-6;   % [F]
R    =  0.010;   % [ohm]

% PI controller parameters (kp + kp/Ti*s)
kpI = 0.427;
kpU = 0.2765;

%%% system equations
% dudt = 1/C*i
% didt = -1/L*u-R/L*i+1/L*u0
% iRef = kpU*(uRef-u)
% u0 = kpI*(iRef-i) = -kpI*i + kpI*(kpU*(uRef-u))

%%% state: x = [u i]
% dxdt = A*x+B*u
%    y = C*x

A = [ 0,            1/Ccap;
    -1/L-kpI*kpU/L, -R/L-kpI/L];
    
 B = [ 0;           kpI*kpU/L];
 
 C = [ 1,           0];
 
 sys = ss(A, B, C, 0);
 
 bode(sys), grid on
 
 pole(sys)
 
 figure
 pzmap(sys), axis equal
end
