
%close all
clear

% figure options
myFontSize = 8;
myLineWidth = 2;
figSize = [530    55   360   380];

FACT = 1;

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
 VSI.tk = VSI.tk -0.0005 - 0.0004 + 0.0015;
 VSI.ton   = 0.1;
% overwrite connection time of load
%VSI.tk    =   10;     % [s] disconnection of the UPS's loads 
%VSI.ton  =    0.3;      % [s]  connection of the UPS's loads 
    
 
VSI.Treg = 30e-06;
%UPS.LoadR123=   1;                                     
%UPS.LoadD3  =   0;                                                  
%UPS.LoadD1  =   0;
VSI.DR_Vc0    = 280;
        
 
%lam1_pr = 111.1867;
%lam3_pr =  313.8555;
%lam5_pr = 344.3558;
%lam7_pr = 690.0204;
%lam9_pr = 0;


lam1_pr = 263.0489;
lam3_pr = 600.0000;
lam5_pr = 294.1350;
lam7_pr = 969.2109;
lam9_pr = 1e-3;

% opt. with new LQR gain and two PRs
lam1_pr = 328.9629;
lam3_pr = 781.0001;
lam5_pr = 0;
lam7_pr = 0;
lam9_pr = 0;
% opt. with new LQR gain and two PRs
lam1_pr = 201.1867;
lam3_pr = 409.8555;
lam5_pr = 346.3558;
lam7_pr = 640.0204;
 
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

    CTRL.s5 = 0;
    
    simOut = sim('UPS_3level_Ttype_Inverter_simplified_ctrl_v3');
        
    load vCf
    load iLf
    load iLoad
    load v_ref
    load i_ref
    
   
Tg = 0.02;
t1 = Tfinal - Tg; 
t2 = Tfinal; 

Imax = 434.7826;



%plot(vCf.Time,v_ref.Data(:,1)/(S.Vn*sqrt(2)),'k'),hold on
%plot(vCf.Time,vCf.Data(:,1)/(S.Vn*sqrt(2)),'b')
%break

% output voltage
t_vec = vCf.Time;
vCf_vec_a = vCf.Data(:,1)/(S.Vn*sqrt(2)*FACT); 
vCf_vec_b = vCf.Data(:,2)/(S.Vn*sqrt(2)*FACT); 
vCf_vec_c = vCf.Data(:,3)/(S.Vn*sqrt(2)*FACT); 
vCf_vec_ref_a = v_ref.Data(:,1)/(S.Vn*sqrt(2)*FACT); 
vCf_vec_ref_b = v_ref.Data(:,2)/(S.Vn*sqrt(2)*FACT); 
vCf_vec_ref_c = v_ref.Data(:,3)/(S.Vn*sqrt(2)*FACT); 

[val1,ind1] = min(abs(t_vec - t1));
[val2,ind2] = min(abs(t_vec - 2));
t_vec_per = t_vec(ind1:ind2)-t_vec(ind1);
vCf_vec_per_a = vCf_vec_a(ind1:ind2,1); 
vCf_vec_per_b = vCf_vec_b(ind1:ind2,1); 
vCf_vec_per_c = vCf_vec_c(ind1:ind2,1); 
vCf_vec_ref_per_a = vCf_vec_ref_a(ind1:ind2,1); 
vCf_vec_ref_per_b = vCf_vec_ref_b(ind1:ind2,1); 
vCf_vec_ref_per_c = vCf_vec_ref_c(ind1:ind2,1); 

h11 = plot(t_vec,vCf_vec_a,'b'),grid on,hold on
plot(t_vec,vCf_vec_b,'b')
plot(t_vec,vCf_vec_c,'b')
plot(t_vec,vCf_vec_ref_a,'k--')
plot(t_vec,vCf_vec_ref_b,'k--')
plot(t_vec,vCf_vec_ref_c,'k--')



lam1_pr = 1e-3;
lam3_pr = 1e-3;
lam5_pr = 1e-3;
lam7_pr = 1e-3;
lam9_pr = 1e-3;



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

    CTRL.s5 = 0;
    
    simOut = sim('UPS_3level_Ttype_Inverter_simplified_ctrl_v3');
        
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



h12 = plot(t_vec,vCf_vec_a,'r'),grid on,hold on
plot(t_vec,vCf_vec_b,'r')
plot(t_vec,vCf_vec_c,'r')
plot(t_vec,vCf_vec_ref_a,'k--')
plot(t_vec,vCf_vec_ref_b,'k--')
plot(t_vec,vCf_vec_ref_c,'k--')
del = 0.007;
%axis([0.4059-del,0.4059+del,-1.5,1.5])
axis([0.4059-del,0.46,-1.5,1.5])

xlabel('time [s]'),ylabel('voltage [pu]')

set(gca,'FontSize', myFontSize);
hleg1 = legend([h11,h12],{'PR','no PR'},'Location','SouthEast')
set(hleg1,'Interpreter','latex')
set(gca,'fontsize',6)
 
break

set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)


matlabfrag('transient_vCf')
movefile('transient_vCf.*', '../fig', 'f')

