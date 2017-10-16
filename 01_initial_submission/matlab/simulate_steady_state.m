
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
% n_T_PWM_1 = 8378; % number of PWM sampling periods when switch in load
% VSI.tk    = n_T_PWM_1*VSI.Tinp/2 + 1e-6 - VSI.Treg - 0.02*30; % [s] disconnection of the UPS's loads 
% VSI.ton   = 0.1;
% overwrite connection time of load
VSI.tk    =   10;     % [s] disconnection of the UPS's loads 
VSI.ton  =    0.11;      % [s]  connection of the UPS's loads 
    
 
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
vCf_vec_ref_a = v_ref.Data(:,1)/(S.Vn*sqrt(2)); 
vCf_vec_ref_b = v_ref.Data(:,2)/(S.Vn*sqrt(2)); 
vCf_vec_ref_c = v_ref.Data(:,3)/(S.Vn*sqrt(2)); 

[val1,ind1] = min(abs(t_vec - t1));
[val2,ind2] = min(abs(t_vec - 2));
t_vec_per = t_vec(ind1:ind2)-t_vec(ind1);
vCf_vec_per_a = vCf_vec_a(ind1:ind2,1); 
vCf_vec_per_b = vCf_vec_b(ind1:ind2,1); 
vCf_vec_per_c = vCf_vec_c(ind1:ind2,1); 
vCf_vec_ref_per_a = vCf_vec_ref_a(ind1:ind2,1); 
vCf_vec_ref_per_b = vCf_vec_ref_b(ind1:ind2,1); 
vCf_vec_ref_per_c = vCf_vec_ref_c(ind1:ind2,1); 

% plot steady state
figure(1)
set(gcf,'outerposition', figSize, 'PaperPositionMode', 'auto')



plot(t_vec_per,vCf_vec_per_a,'b'),grid on,hold on
plot(t_vec_per,vCf_vec_per_b,'r')
plot(t_vec_per,vCf_vec_per_c,'g')
plot(t_vec_per,vCf_vec_ref_per_a,'k--')
plot(t_vec_per,vCf_vec_ref_per_b,'k--')
plot(t_vec_per,vCf_vec_ref_per_c,'k--')
legend('a','b','c','Location','South', 'orientation', 'horizontal')
xlabel('time [s]'),ylabel('phase voltage [V]'), ylim([-1.4, 1.1])
set(gca,'FontSize', myFontSize);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)
matlabfrag('steady_state_vCf')

movefile('steady_state_vCf.*', '../fig', 'f')
 

% load current
iLoad_vec_a = iLoad.Data(:,1)/Imax; 
iLoad_vec_b = iLoad.Data(:,2)/Imax; 
iLoad_vec_c = iLoad.Data(:,3)/Imax; 

iLoad_vec_per_a = iLoad_vec_a(ind1:ind2,1); 
iLoad_vec_per_b = iLoad_vec_b(ind1:ind2,1); 
iLoad_vec_per_c = iLoad_vec_c(ind1:ind2,1); 


figure(3)
set(gcf,'outerposition', figSize, 'PaperPositionMode', 'auto')

plot(t_vec_per,iLoad_vec_per_a,'b'),grid on,hold on
plot(t_vec_per,iLoad_vec_per_b,'r')
plot(t_vec_per,iLoad_vec_per_c,'g')
legend('a','b','c','Location','South', 'orientation', 'horizontal')
ylabel('phase current [A]'), ylim([-1.8, 1.4])
xlabel('time [s]'), set(gca,'FontSize', myFontSize);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)
matlabfrag('steady_state_iLoad')

movefile('steady_state_iLoad.*', '../fig', 'f')
 



vCf_FFT_a = fft(vCf_vec_per_a)/length(vCf_vec_per_a);
vCf_THD_a = 100 * sqrt(sum((2*abs(vCf_FFT_a(4:2:10))).^2))


% figure(4)
% set(gcf,'outerposition', figSize, 'PaperPositionMode', 'auto')
% 
% bar(20*log10(2*abs(vCf_FFT_a(1:50)))),grid on
% axis([0,50,-60,0])
% 
% xlabel('harmonic number'),ylabel('amplitude [dB]')
% set(gca,'FontSize', myFontSize);
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)
% matlabfrag('steady_state_harmonics')
% 
% movefile('steady_state_harmonics.*', '../fig', 'f')
 
figure(5)
set(gcf,'outerposition', figSize, 'PaperPositionMode', 'auto')

bar(0:49, 2*abs(vCf_FFT_a(1:50))),grid on
axis([-1,50,0,0.05])
text(2,0.048, '$h_1 = 0.9946$')

xlabel('harmonic number'),ylabel('amplitude [pu]')
set(gca,'FontSize', myFontSize);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)
matlabfrag('steady_state_harmonics')

movefile('steady_state_harmonics.*', '../fig', 'f')








