
close all
clear

% figure options
myFontSize = 8;
myLineWidth = 2;
figSize = [530    55   360   380];

% ------------------------------------------------------------- 
% The controller varieties considered in simulation
% are summarised in the report

Run

SW_POS = [1 0 0 0 0 0;
          0 1 1 0 0 0;
          0 0 1 1 0 0;
          0 1 1 0 1 0;
          0 1 1 0 1 1;
          0 0 1 1 1 0;
          0 0 1 1 1 1];

      % test cases:
% load_conf 1; single resisor connected to A2 phase
% load_conf 2; two resistors connected to A2 and B2 phase, respectively                                                                   
% load_conf 3  three resistors (balansed)
% load_conf 4  three phase diode rectifier
% load_conf 5  single phase diode rectifier

load_conf = 4;


% overwrite type of load
switch load_conf
    case 1
       UPS.LoadR123=   1;                                       
       UPS.LoadD3  =   0;                                                  
       UPS.LoadD1  =   0;     
    case 2
       UPS.LoadR123=   2;                                            
       UPS.LoadD3  =   0;                                                 
       UPS.LoadD1  =   0;     
    case 3
       UPS.LoadR123=   3;                                         
       UPS.LoadD3  =   0;                                                  
       UPS.LoadD1  =   0;     
    case 4
       UPS.LoadR123=   0;                                         
       UPS.LoadD3  =   1;                                                  
       UPS.LoadD1  =   0;     
    case 5
       UPS.LoadR123=   0;                                         
       UPS.LoadD3  =   0;                                                  
       UPS.LoadD1  =   1;        
end


 Tfinal = 0.5; % 1.1;
 n_T_PWM_1 = 8378; % number of PWM sampling periods when switch in load
 VSI.tk    = n_T_PWM_1*VSI.Tinp/2 + 1e-6 - VSI.Treg - 0.02*30; % [s] disconnection of the UPS's loads 
 VSI.ton   = 0.1;
 VSI.tk    = VSI.tk - 0.0005;
 
VSI.Treg = 30e-06;
%UPS.LoadR123=   1;                                     
%UPS.LoadD3  =   0;                                                  
%UPS.LoadD1  =   0;
VSI.DR_Vc0    = 280;
        
LQR_control_design


    % --- run simulations -----------------------------------
    % -------------------------------------------------------
    ctrl_case = 7 %1,...,7

    CTRL.s0 = SW_POS(ctrl_case,1);
    CTRL.s1 = SW_POS(ctrl_case,2);
    CTRL.s2 = SW_POS(ctrl_case,3);
    CTRL.s3 = SW_POS(ctrl_case,4);
    CTRL.s4 = SW_POS(ctrl_case,5);
    CTRL.s5 = SW_POS(ctrl_case,6);

    simOut = sim('UPS_3level_Ttype_Inverter_simplified_ctrl');
        
    load vCf
    load iLf
    load iLoad
    load v_ref
    load i_ref
    
   
Tg = 0.02;
t1 = Tfinal - Tg; 
t2 = Tfinal; 

Imax = 434.7826;





% output voltage
t_vec = vCf.Time;
vCf_vec_a = vCf.Data(:,1)/(S.Vn*sqrt(2)); 
vCf_vec_b = vCf.Data(:,2)/(S.Vn*sqrt(2)); 
vCf_vec_c = vCf.Data(:,3)/(S.Vn*sqrt(2)); 


figure(1)
set(gcf,'outerposition', figSize, 'PaperPositionMode', 'auto')

plot(t_vec, vCf_vec_a,'b'),grid on,,hold on
plot(t_vec, vCf_vec_b,'r')
plot(t_vec, vCf_vec_c,'g')
axis([0.395,0.415,-1.5,1.5])

legend('a','b','c','Location','South', 'orientation', 'horizontal')
xlabel('time [s]'),ylabel('phase voltage [V]'), ylim([-1.5, 1.3])
set(gca,'FontSize', myFontSize);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)
matlabfrag('transient_vCf')

movefile('transient_vCf.*', '../fig', 'f')

