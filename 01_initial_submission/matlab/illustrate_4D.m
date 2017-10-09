

%yalmip('clear')
%oldopts = sdpsettings;
%solveropts = sdpsettings(oldopts,'verbose',0,'solver','sedumi');


options = optimset;
options = optimset(options,'Display','off');


L =  85e-6;   % [H]
Ccap = 275e-6;   % [F]
R = 0.010;    % [ohm]

% G = b/(s^2 a*s + b)
b = 10463/(Ccap*10000*L);
a = (3259 + 10000*R)/(10000*L);
w1 = 2*pi*50;
w3 = 3*w1;
w5 = 5*w1;
w7 = 7*w1;

opt_gains = [  111.1867;
  313.8555;
  344.3558;
  690.0204];
  
% -----------------------------------------------------
% least squares fitting over grid
% lower, upper bound of gains
L = 0;
U = 1000;
C = [-1; 0; 0; 0; 0];
% grid size of boxes 
h1 = 100; % 100; % box size
% enumerate the boxes
LV = L:h1:U-h1; 
figure(1),hold on
for kkk1 = 1:length(LV)
for kkk2 = 1:length(LV)   
for kkk3 = 1:length(LV)
     [kkk1, kkk2, kkk3]
for kkk4 = 1:length(LV)
    
        L1 = LV(kkk1);
        L3 = LV(kkk2); 
        L5 = LV(kkk3); 
        L7 = LV(kkk4); 
        xk = [L1; L3; L5; L7];
        
        den = [ 1;
 a ;
 (w1^2 + w3^2 + w5^2 + w7^2 + b);
 (a*w1^2 + a*w3^2 + a*w5^2 + a*w7^2 + L1*b + L3*b + L5*b + L7*b) ;
 (w1^2*w3^2 + w1^2*w5^2 + w1^2*w7^2 + b*w1^2 + w3^2*w5^2 + w3^2*w7^2 + b*w3^2 + w5^2*w7^2 + b*w5^2 + b*w7^2);
 (a*w1^2*w3^2 + a*w1^2*w5^2 + a*w1^2*w7^2 + a*w3^2*w5^2 + a*w3^2*w7^2 + a*w5^2*w7^2 + L1*b*w3^2 + L3*b*w1^2 + L1*b*w5^2 + L5*b*w1^2 + L1*b*w7^2 + L3*b*w5^2 + L5*b*w3^2 + L7*b*w1^2 + L3*b*w7^2 + L7*b*w3^2 + L5*b*w7^2 + L7*b*w5^2);
 (w1^2*w3^2*w5^2 + w1^2*w3^2*w7^2 + b*w1^2*w3^2 + w1^2*w5^2*w7^2 + b*w1^2*w5^2 + b*w1^2*w7^2 + w3^2*w5^2*w7^2 + b*w3^2*w5^2 + b*w3^2*w7^2 + b*w5^2*w7^2);
 (a*w1^2*w3^2*w5^2 + a*w1^2*w3^2*w7^2 + a*w1^2*w5^2*w7^2 + a*w3^2*w5^2*w7^2 + L1*b*w3^2*w5^2 + L3*b*w1^2*w5^2 + L5*b*w1^2*w3^2 + L1*b*w3^2*w7^2 + L3*b*w1^2*w7^2 + L7*b*w1^2*w3^2 + L1*b*w5^2*w7^2 + L5*b*w1^2*w7^2 + L7*b*w1^2*w5^2 + L3*b*w5^2*w7^2 + L5*b*w3^2*w7^2 + L7*b*w3^2*w5^2);
 (w1^2*w3^2*w5^2*w7^2 + b*w1^2*w3^2*w5^2 + b*w1^2*w3^2*w7^2 + b*w1^2*w5^2*w7^2 + b*w3^2*w5^2*w7^2) ;
 (a*w1^2*w3^2*w5^2*w7^2 + L7*b*w1^2*w3^2*w5^2 + L5*b*w1^2*w3^2*w7^2 + L3*b*w1^2*w5^2*w7^2 + L1*b*w3^2*w5^2*w7^2);
 b*w1^2*w3^2*w5^2*w7^2  ];

        rot = roots(den);
        plot(rot,'g.')
        
end
end
end
end


L1 = opt_gains(1);
L3 = opt_gains(2); 
L5 = opt_gains(3); 
L7 = opt_gains(4); 
       
        
        den = [ 1;
 a ;
 (w1^2 + w3^2 + w5^2 + w7^2 + b);
 (a*w1^2 + a*w3^2 + a*w5^2 + a*w7^2 + L1*b + L3*b + L5*b + L7*b) ;
 (w1^2*w3^2 + w1^2*w5^2 + w1^2*w7^2 + b*w1^2 + w3^2*w5^2 + w3^2*w7^2 + b*w3^2 + w5^2*w7^2 + b*w5^2 + b*w7^2);
 (a*w1^2*w3^2 + a*w1^2*w5^2 + a*w1^2*w7^2 + a*w3^2*w5^2 + a*w3^2*w7^2 + a*w5^2*w7^2 + L1*b*w3^2 + L3*b*w1^2 + L1*b*w5^2 + L5*b*w1^2 + L1*b*w7^2 + L3*b*w5^2 + L5*b*w3^2 + L7*b*w1^2 + L3*b*w7^2 + L7*b*w3^2 + L5*b*w7^2 + L7*b*w5^2);
 (w1^2*w3^2*w5^2 + w1^2*w3^2*w7^2 + b*w1^2*w3^2 + w1^2*w5^2*w7^2 + b*w1^2*w5^2 + b*w1^2*w7^2 + w3^2*w5^2*w7^2 + b*w3^2*w5^2 + b*w3^2*w7^2 + b*w5^2*w7^2);
 (a*w1^2*w3^2*w5^2 + a*w1^2*w3^2*w7^2 + a*w1^2*w5^2*w7^2 + a*w3^2*w5^2*w7^2 + L1*b*w3^2*w5^2 + L3*b*w1^2*w5^2 + L5*b*w1^2*w3^2 + L1*b*w3^2*w7^2 + L3*b*w1^2*w7^2 + L7*b*w1^2*w3^2 + L1*b*w5^2*w7^2 + L5*b*w1^2*w7^2 + L7*b*w1^2*w5^2 + L3*b*w5^2*w7^2 + L5*b*w3^2*w7^2 + L7*b*w3^2*w5^2);
 (w1^2*w3^2*w5^2*w7^2 + b*w1^2*w3^2*w5^2 + b*w1^2*w3^2*w7^2 + b*w1^2*w5^2*w7^2 + b*w3^2*w5^2*w7^2) ;
 (a*w1^2*w3^2*w5^2*w7^2 + L7*b*w1^2*w3^2*w5^2 + L5*b*w1^2*w3^2*w7^2 + L3*b*w1^2*w5^2*w7^2 + L1*b*w3^2*w5^2*w7^2);
 b*w1^2*w3^2*w5^2*w7^2  ];

        rot = roots(den);
        plot(rot,'ro','MarkerSize',3)

        
        
fontsize = 7;
set( 0, 'defaultAxesFontSize', fontsize);
        
 plot([0,real(rot(1))],[0,imag(rot(1))],'r')
 plot([0,real(rot(2))],[0,imag(rot(2))],'r')
 axis([-2000,100,-8000,8000])
 %figuresize(7.2,5.2,'cm')
 % figuresize(8.2,5.6,'cm')
 print -depsc root_locus
       
 
 
 
        
        
        
