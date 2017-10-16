     
% Run.m file allows to define all parameters which are necessary to simulate 2-level, 3- and 4-leg topologies of UPS system in MATLAB/Simulink.
% The file is divided into nine sections marked as (1),(2),(3),(4),(5),(6),(7),(8) and (9).
% Section (1) relates to 


    %close all;
    %clear all;
    %clc;

    global Ts UPS S VD PR VB VSI ADC        % global variables
    
        
%% (1) 

    Ts          =   1/(1*1e6);      % [s]   fixed-step size (fundamental sample time) in MATLAB/Simulink

    Pmax        =   100e3;          % [W]   nominal (max) power of UPS       
    Vsf         =   230;            % [V]   nominal input and output AC voltage
    fn          =   50;             % [Hz]      
    Vd          =   720;            % [V]       
    amVd        =   15;             % [V]       
    Vdmax       =   Vd + amVd;      % [V]       
    Vdmin       =   Vd - amVd;      % [V]                                      
    Imax        =   3.0*Pmax/(3*Vsf);   % [A]            
     

    Vdr         =   sqrt((Vdmax^2 + Vdmin^2)/2);                    % [V]
    Cd          =   27e-3;%   % [F]   
    

    
    S.Vhon      =   1;      % Simulate harmonics in grid voltage (on/off)
                            % information about magnitudes and phases of higher harmonics is placed in section 2 and S.x structure  
                                                                    
    VD.VDx      =   0;      % Simulate input voltage sags (on/off) 
                            % settings of type of voltage sag, its start as well as duratioan are available in section 3 and VD.x structure                         
                                                                    
    UPS.Lfnonlin   =   1;   % Use non-linear inductors (on/off)
                            % information about the inductors' characteristics is placed in: 
                            % - section 4 and PR.x structure for input inductors (Lf1) 
                            % - section 6 and VSI.x structure for output inductors (Lf2)
                                                                                                                                                                                                     
    UPS.LoadR123=   3;      % 0 - does not allow to simulate unbalanced, balanced and linear load by means of resistors 
                            % 1 - allows to simulate unbalanced and linear load by means of a single resisor connected to A2 phase
                            % 2 - allows to simulate unbalanced and linear load by means of two resistors connected to A2 and B2 phase, respectively                                                                   
                            % 3 - allows to simulate balanced and linear load by means of three resistors
                            % resistors' parameters are available in section 6 and VSI.x structure
                                                                    
    UPS.LoadD3  =   0;      % connect 3-phase diode rectifier (on/off)
                                                                    
    UPS.LoadD1  =   0;      % connect 1-phase diode rectifier load (on/off)

    UPS.Balancer = 1;       % activate balancer (on/off)
                            
                            
%% (2) (S.x) structure System parameters

   S.Vn         =   Vsf;            % [V]    RMS value of phase voltage (vg1a,vg1b,vg1c)
   S.Vnpp       =   S.Vn*sqrt(3);   % [V]    RMS value of phase-to-phase voltage (vg1ab,vg1bc,vg1ca)                          
   S.fn         =   fn;             % [Hz]   inpute frequency
   S.Imax       =   Imax;           % [A]    RMS value of maximum input current (ig1a,ig1b,ig1c)
   

   S.Rs         =   0.005;                              % [ohm]  input resistance (Rg1)
   S.Ls         =   71.5e-6;                              % [H]    input inductance (Lg1)
   S.Rb         =   S.Vn*sqrt(2)/(S.Imax*sqrt(2));      % [ohm]  input resistance (Rb1) limiting the capacitor Cd initial charging current      

                                                             
   S.Vh1        =   1*(310.555)/(sqrt(2)*Vsf);     % [pu]  magnitude of 1st harmonic (h=1) in inpute voltage (vg1a,vg1b,vg1c)         
   S.fi0h1      =   262.5613*pi/180;               % [rad] phase of 1st harmonic (h=1) in inpute voltage (vg1a,vg1b,vg1c)                     
            
   S.Vh3        =   S.Vhon*(3.6)/(sqrt(2)*Vsf);    % [pu]  magnitude of 3rd harmonic (h=3) in inpute voltage (vg1a,vg1b,vg1c)
   S.fi0h3      =   24.0124*pi/180;                % [rad] phase of 3rd harmonic (h=3) in inpute voltage (vg1a,vg1b,vg1c)
   
   S.Vh5        =   S.Vhon*(6.15)/(sqrt(2)*Vsf);   % [pu]  magnitude of 5th harmonic (h=5) in inpute voltage (vg1a,vg1b,vg1c)
   S.fi0h5      =   36.6192*pi/180;                % [rad] phase of 5th harmonic (h=5) in inpute voltage (vg1a,vg1b,vg1c)
   
   S.Vh7        =   S.Vhon*(4.78)/(sqrt(2)*Vsf);   % [pu]  magnitude of 7th harmonic (h=7) in inpute voltage (vg1a,vg1b,vg1c)         
   S.fi0h7      =   241.8570*pi/180;               % [rad] phase of 7th harmonic (h=7) in inpute voltage (vg1a,vg1b,vg1c)
   
   S.Vh9        =   S.Vhon*(1.2)/(sqrt(2)*Vsf);  % [pu]  magnitude of 9th harmonic (h=9) in inpute voltage (vg1a,vg1b,vg1c)
   S.fi0h9      =   31.0122*pi/180;                % [rad] phase of 9th harmonic (h=9) in inpute voltage (vg1a,vg1b,vg1c)
   
   S.Vh11       =   S.Vhon*(0.8)/(sqrt(2)*Vsf);   % [pu]  magnitude of 11th harmonic (h=11) in inpute voltage (vg1a,vg1b,vg1c)
   S.fi0h11     =   72.3976*pi/180;                % [rad] phase of 11th harmonic (h=11) in inpute voltage (vg1a,vg1b,vg1c)
   
   S.Vh13       =   S.Vhon*(0.9)/(sqrt(2)*Vsf); % [pu]  magnitude of 13th harmonic (h=13) in inpute voltage (vg1a,vg1b,vg1c)
   S.fi0h13     =   193.2643*pi/180;               % [rad] phase of 13th harmonic (h=13) in inpute voltage (vg1a,vg1b,vg1c)
            
   S.THDvg1     =   (sqrt(S.Vh3^2 + S.Vh5^2 + S.Vh7^2 + S.Vh9^2 + S.Vh11^2 + S.Vh13^2)/S.Vh1)*100;   % [%] THD factor for the input voltage

   
   S.t0         =   0.03;       % [s] UPS starts work - Cd capacitor is charged by Rb resistor   
   S.tk         =   0.10;       % [s] Rb resistor is bypassed

%% (3) (VD.x) Structure, Voltage sag

    VD.t0       =   0.45;     % [s] start of inpute voltage sag (vg1a,vg1b,vg1c) to a value fixed in the VD.VD1 table (see below) 
                              %     provided that the VD.VDx coefficient equals 1    
    VD.tVD      =   0.1;      % [s] duration of the voltage sag
    VD.VD1      =   [1    0.40    0.40];     % [pu] amplitude of the input voltage (vg1a,vg1b,vg1c) when the voltage sag is set (1 = S.Vh1) 

%% (4) (PR.x) Structure, Rectifier settings

    L2NonLin_I_PSI=load('L2.1_N30_3co__I_Psi.txt');

    PR.Vdr          =   Vdr;          % [V] - reference value of the DC-link voltage (Vdc), (see section 1)
    PR.t0           =   0.05;         % [s] rectifier starts boosting the Vdc voltage to the reference value 
  
    PR.DC_Rload     =   PR.Vdr^2/Pmax; % dc-load in case of rectifier only simulation
    
    PR.finp         =   4200;         % [Hz] - switching frequency of the rectifier's PWM modulator (fsw)
    PR.Tinp         =   1/PR.finp;
       
    PR.freg         =   33600;        % [Hz] - sample frequency of the rectifier's control system. In Simulink, the rectifier's control system   
    PR.Treg         =   1/PR.freg;    %        is placed in a block signed as 'PFC Control System SRF(dq0) - SPWM' 
   
    PR.Cd           =   Cd            % [F]   - capacity of the DC-link capacitor (Cdc)
    PR.Rcd          =   0.01;         % [ohm] - series resisance of the DC-link capacitor (Cdc) 
    PR.Vcd0         =   650;          % [V]   - initial voltage of the DC-link capacitor (Cdc)
 
    PR.Rcf          =   0.005;        % [ohm] - series resistance of the input filter's capacitors (Cf1)
    PR.Cf           =   275e-6;       % [F]   - capacity of the input filter's capacitors (Cf1)
    PR.Vcf0         =   0;            % [V]   - initial voltage of the input filter's capacitors (Cf1)
    
    
    PR.Rf           =   0.050;                                              % [ohm] - series resistance of the inductor Lf1 
    PR.RfFe         =   100;                                                % [ohm] - parallel resistance of the inductor Lf1
    PR.Lf           =   94.8e-6;                                            % [H]   - inductance of the inductor Lf1
    PR.IF1          =   L2NonLin_I_PSI(1,:);                                %       - necessary vectors for simulation of the linear inductance Lf1 (=171e-6 [H])     
    PR.FLUXlin1     =   PR.Lf*PR.IF1;                                       %         by means of Lookup Table block (see the input filter in Simulink)
    PR.IF2          =   -fliplr(PR.IF1(2:end));                             %               - || -
    PR.FLUXlin2     =   -fliplr(PR.FLUXlin1(2:end)) + 2*PR.FLUXlin1(1);     %               - || -
    PR.IF           =   [PR.IF2 PR.IF1];                                    %               - || -
    PR.FLUXlin      =   [PR.FLUXlin2  PR.FLUXlin1];                         %               - || -
    
    
    PR.FLUXnonlin1  =   L2NonLin_I_PSI(2,:);
    PR.FLUXnonlin2  =   -fliplr(PR.FLUXnonlin1(2:end)) + 2*PR.FLUXnonlin1(1);   %    by means of Lookup Table block (see the input filter in Simulink)
    PR.FLUXnonlin   =   [PR.FLUXnonlin2  PR.FLUXnonlin1];                       %           - || -
    PR.Lfnonlin     =   PR.FLUXnonlin1(100:end)./PR.IF1(100:end);               %           - || -
    
    
% dc-link voltage control
    PR.Kp1          =   6.25;               %     - proportional constant of the rectifier's voltage controller (external controller) 
    PR.Ti1          =   0.025;              %     - integral constant of the rectifier's voltage controller 
    PR.Tsam1        =   PR.Treg;            % [s] - sample time of the rectifier's voltage controller 
    PR.lim_up1      =   1*sqrt(2)*Imax;     % [A] - upper limitation of the rectifier's voltage controller 
    PR.lim_low1     =   -1*sqrt(2)*Imax;    % [A] - lower limitation of the rectifier's voltage controller 

% input current control
    PR.Kp2       =   0.4273;                %     - proportional constant of the rectifier's current controllers (internal controllers)  
    PR.Ti2       =   0.0047;                %     - integral constant of the rectifier's current controllers 
    PR.Tsam2     =   PR.Treg;               % [s] - sample time of the rectifier's current controllers
    PR.lim_up2   =   1*sqrt(2)*S.Vn;        % [V] - upper limitation of the rectifier's current controllers
    PR.lim_low2  =   -1*sqrt(2)*S.Vn;       % [V] - lower limitation of the rectifier's current controllers
    
  
    PR.ttransf   =   1*1e-6;                % [s] - delay generated by drivers of the output IGBT transistors 
    PR.Nttransf  =   floor(PR.ttransf/Ts);  %           - || -
    
    PR.tdead     =   3*1e-6;                % [s] - dead time of the input IGBT transistors 
    PR.Ntdead    =   floor(PR.tdead/Ts);    %           - || -
    
    
       
%% (5) (VB.x) Structure
        
    VB.finp      =  10000;              % [Hz] - switching frequency of the voltage balancing circuit (VBC)
    VB.Tinp      =  1/VB.finp;      
    

    VB.Kp1       =   5;                 %     - proportional constant of the VBC's voltage controller (external controller) 
    VB.Ti1       =   0.1;               %     - integral constant of the VBC's voltage controller 
    VB.Tsam1     =   VB.Tinp/2;         % [s] - sample time of the VBC's voltage controller 
    VB.lim_up1   =   1*Imax*sqrt(2);    % [A] - upper limitation of the VBC's voltage controller 
    VB.lim_low1  =   -1*Imax*sqrt(2);   % [A] - lower limitation of the VBC's voltage controller
    
    
    VB.Kp2       =   0.1;               %     - proportional constant of the VBC's current controller (internal controller)             
    VB.Ti2       =   0.02;              %     - integral constant of the VBC's current controller
    VB.Tsam2     =   VB.Tinp;           % [s] - sample time of the VBC's current controller 
    VB.lim_up2   =   360;               %     - upper limitation of the VBC's current controller 
    VB.lim_low2  =   -360;              %     - lower limitation of the VBC's current controller
    
    
    VB.R         =   0.05;      % [ohm] - series resistance of the VBC's inductors (Lvb1,Lvb2)
    VB.L         =   200e-6;    % [H]   - inductace of the VBC's inductors (Lvb1,Lvb2) 
      
    VB.ttransf   =   1*1e-6;                % [s] - delay generated by drivers of the IGBT transistors
    VB.Nttransf  =   floor(VB.ttransf/Ts);  %           - || -
    
    VB.t0        =   0;                 % starting time of voltage balancer
    
%% (6) (VSI.x) Structure

    VSI.finp        =   4200;           % [Hz] - switching frequency of the output inverter's PWM modulator (fsw) 
    %--%VSI.Tinp        =   1/VSI.finp;
    % make switching period an integer 
    % multiple of the simulating tim 1e-6 s
    %--%VSI.Tinp        = 238e-06;
    VSI.Tinp        = 240e-06; % take this to make it a multiple of control period
    
    
    VSI.freg        =   33600;          % [Hz] - sample frequency of the output inverter's control system. In Simulink, the inverter's control system
    %--%VSI.Treg        =   1/VSI.freg;     %        is placed in a block signed as 'Output Inverter Control System SRF(dq0) - SPWM' 
    % make switching period an integer 
    % multiple of the simulating tim 1e-6 s
    VSI.Treg        = 30e-06;
    
    VSI.kVd         =   8;
    
    VSI.Rcf         =   0.005;          % [ohm] - series resistance of the output filter's capacitors (Cf2)
    VSI.Cf          =   275e-6;         % [F]   - capacity of the output filter's capacitors (Cf2)
    VSI.Vcf0        =   0;              % [V]   - initial voltage of the output filter's capacitors (Cf2)
    VSI.Icf         =   32.5;
    
   
    VSI.Rf          =   0.010;                                              % [ohm] - series resistance of the inductor Lf2 
    VSI.RfFe        =   1000;                                               % [ohm] - parallel resistance of the inductor Lf2
    VSI.Lf          =   94.8e-6;                                            % [H]   - inductance of the inductor Lf2
    VSI.IF1         =   L2NonLin_I_PSI(1,:);                                %       - necessary vectors for simulation of the linear inductance Lf2 (=171e-6 [H]) 
    VSI.FLUXlin1    =   VSI.Lf*VSI.IF1;                                     %         by means of Lookup Table block (see the output filter in Simulink)
    VSI.IF2         =   -fliplr(VSI.IF1(2:end));                            %               - || -
    VSI.FLUXlin2    =   -fliplr(VSI.FLUXlin1(2:end)) + 2*VSI.FLUXlin1(1);   %               - || -
    VSI.IF          =   [VSI.IF2 VSI.IF1];                                  %               - || -
    VSI.FLUXlin     =   [VSI.FLUXlin2  VSI.FLUXlin1];                       %               - || -
    VSI.Lflin       =   VSI.FLUXlin1(2:end)./VSI.IF1(2:end);                %               - || -
    
    
    VSI.FLUXnonlin1 =   L2NonLin_I_PSI(2,:);                                        % -  necessary vectors for simulation of the nonlinear inductance Lf2
    VSI.FLUXnonlin2 =   -fliplr(VSI.FLUXnonlin1(2:end)) + 2*VSI.FLUXnonlin1(1);     %    by means of Lookup Table block (see the output filter in Simulink)
    VSI.FLUXnonlin  =   [VSI.FLUXnonlin2  VSI.FLUXnonlin1];                         %           - || -
    VSI.Lfnonlin    =   VSI.FLUXnonlin1(100:end)./VSI.IF1(100:end);                 %           - || -

    VSI.fres      =   1/(2*pi*sqrt(VSI.Lf*VSI.Cf));        % [Hz] - resonant frequency of the output filter (Cf2,Lf2)
    
% output voltage controller
    VSI.Kp1       =   0.2765;                              %     - proportional constant of the inverter's voltage controllers (external controllers) 
    VSI.Ti1       =   0.036;                               %     - integral constant of the inverter's voltage controllers 
    VSI.Tsam1     =   VSI.Treg;                            % [s] - sample time of the inverter's voltage controllers 
    VSI.lim_up1   =   1*sqrt(2)*sqrt(VSI.Icf^2 + Imax^2)   % [A] - upper limitation of the inverter's voltage controllers
    VSI.lim_low1  =   -1*sqrt(2)*sqrt(VSI.Icf^2 + Imax^2)  % [A] - lower limitation of the inverter's voltage controllers
    
% output current controller
    VSI.Kp2       =   0.4273;           %     - proportional constant of the inverter's current controllers (internal controllers)
    VSI.Ti2       =   0.000665;         %     - integral constant of the inverter's current controllers
    VSI.Tsam2     =   VSI.Treg;         % [s] - sample time of the inverter's current controllers
    VSI.lim_up2   =   1*S.Vn*sqrt(2);   % [V] - upper limitation of the inverter's current controllers
    VSI.lim_low2  =   -1*S.Vn*sqrt(2);  % [V] - lower limitation of the inverter's current controllers

    
    VSI.Ra        = 3*Vsf^2/(1*Pmax);     % [ohm] - resistance of the resistive load for A2 phase
    VSI.La        = 0;                    % [H]   -
    
      
    VSI.Rb        = 3*Vsf^2/(1*Pmax);     % [ohm] - resistance of the resistive load for B2 phase     
    VSI.Lb        = 0;                    % [H]   -
  
   
    VSI.Rc        = 3*Vsf^2/(1*Pmax);     % [ohm] - resistance of the resistive load for C2 phase
    VSI.Lc        = 0;                    % [H]   -
    
          
    VSI.DR_Rs     = 64.1e-3;     % [Ohm] - inductance of the 1- and 3-phase diode rectifier which simulate an unbalanced (balanced) and nonlinera load
    VSI.DR_C      = 41.5e-3;      % [F] - capacity of the 1- and 3-phase diode rectifier
    VSI.DR_Vc0    = 280;          % [V] - initial voltage of the the diode recifiers' capacitors
    VSI.DR_Rload  = 3.615;
    
    VSI.ttransf   =   1*1e-6;                   % [s] - delay generated by drivers of the output IGBT transistors
    VSI.Nttransf  =   floor(VSI.ttransf/Ts);    %           - || -
    
    
    VSI.tdead     =   2*1e-6;                   % [s] - dead time of the output IGBT transistors 
    VSI.Ntdead    =   floor(VSI.tdead/Ts);      %           - || -
   
    
    VSI.t0        =   0.05;%0.18;     % [s] output inverter's algorithm starts work
  
    VSI.tk        =   0.5;            % [s] disconnection of the UPS's loads 
 
    %VSI.ton       =   VSI.tk + 0.2;   % [s] connection of the UPS's loads 
    VSI.ton       =   0.2;   % [s] connection of the UPS's loads 
    
    %% Harmonic compensators
    % Rectifier
    s=tf('s');

    w0HC3=2*pi*150;
    phi3=0.3;
    Tn3=0.2;
    GrHC3=2*1/Tn3*(s*cos(phi3)-w0HC3*sin(phi3))/(s^2+w0HC3^2);

    w0HC5=2*pi*250;
    phi5=0.5;
    Tn5=0.2;
    GrHC5=2*1/Tn5*(s*cos(phi5)-w0HC5*sin(phi5))/(s^2+w0HC5^2);

    w0HC7=2*pi*350;
    phi7=0.7;
    Tn7=0.2;
    GrHC7=2*1/Tn7*(s*cos(phi7)-w0HC7*sin(phi7))/(s^2+w0HC7^2);

    HC3=c2d(GrHC3,PR.Treg,'tustin')
    [numHC3rect,denHC3rect]=tfdata(HC3,'v')
    HC5=c2d(GrHC5,PR.Treg,'tustin')
    [numHC5rect,denHC5rect]=tfdata(HC5,'v')
    HC7=c2d(GrHC7,PR.Treg,'tustin')
    [numHC7rect,denHC7rect]=tfdata(HC7,'v')
  
    
    % Inverter
    w0HC3=2*pi*150;
    phi3=0.3;
    Tn3=0.2;
    GrHC3=2*1/Tn3*(s*cos(phi3)-w0HC3*sin(phi3))/(s^2+w0HC3^2);

    w0HC5=2*pi*250;
    phi5=0.5;
    Tn5=0.2;
    GrHC5=2*1/Tn5*(s*cos(phi5)-w0HC5*sin(phi5))/(s^2+w0HC5^2);

    w0HC7=2*pi*350;
    phi7=0.7;
    Tn7=0.2;
    GrHC7=2*1/Tn7*(s*cos(phi7)-w0HC7*sin(phi7))/(s^2+w0HC7^2);

    HC3=c2d(GrHC3,VSI.Treg,'tustin')
    [numHC3,denHC3]=tfdata(HC3,'v')
    HC5=c2d(GrHC5,VSI.Treg,'tustin')
    [numHC5,denHC5]=tfdata(HC5,'v')
    HC7=c2d(GrHC7,VSI.Treg,'tustin')
    [numHC7,denHC7]=tfdata(HC7,'v')
    

    

    
%% (7) (ADC.x) Structure

    ADC.N         = 12;     %     - number of ADC bits
    ADC.Vmin      = 0;      % [V] - minimum value of an input voltage 
    ADC.Vmax      = 3;      % [V] - maximum value of an input voltage
    ADC.Voffset   = 1.5;    % [V] - offset
    ADC.Nonlin    = 1;      %     - level 
    ADC.Nnoise    = 1;

    
    ADC.k1vdc     = 1.5/1000;           %[V]/[V] - measuring channel's gain
    ADC.k2vdc     = 1000/(2^ADC.N/2);   
    
    ADC.k1vDCpn     = 1.5/500;           %[V]/[V] - measuring channel's gain
    ADC.k2vDCpn     = 500/(2^ADC.N/2); 
    

    ADC.k1vABC     = 1.5/600;           %[V]/[V] - measuring channel's gain
    ADC.k2vABC     = 600/(2^ADC.N/2);
    
    
    ADC.k1iABC     = 1.5/600;           %[V]/[A] - measuring channel's gain 
    ADC.k2iABC     = 600/(2^ADC.N/2);
    
    
    ADC.k1iLb      = 1.5/600;           %[V]/[A] - measuring channel's gain 
    ADC.k2iLb      = 600/(2^ADC.N/2);

%% (8)
 	                                                                
    %FiguresLf1Lf2;
    
%% (9) Update nonlinear current flux relations according to CECUA doc
load NonlinearFluxCurrentRelation

    % load side inductance Lfi
VSI.Lf         =   127.8e-6;
VSI.IF         =   LCL.I;
VSI.FLUXlin    =   VSI.Lf.*VSI.IF;
VSI.FLUXnonlin =   LCL.FLUX;
VSI = rmfield(VSI, {'IF1', 'IF2', 'FLUXlin1', 'FLUXlin2', ...
    'FLUXnonlin1', 'FLUXnonlin2', 'Lfnonlin'});

    
    