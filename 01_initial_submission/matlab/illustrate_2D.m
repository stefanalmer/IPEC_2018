

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
%w5 = 5*w1;
%w7 = 7*w1;

%opt_gains = [  111.1867;
%  313.8555;
%  344.3558;
%  690.0204];
  
% -----------------------------------------------------
% least squares fitting over grid
% lower, upper bound of gains
L = 0;
U = 1000;
C = [-1; 0; 0; 0; 0];
% grid size of boxes 
h1 = 10; % 100; % box size
% enumerate the boxes
LV = L:h1:U-h1; 
figure(1),hold on
for kkk1 = 1:length(LV)
    kkk1
for kkk2 = 1:length(LV)   
%for kkk3 = 1:length(LV)
%     [kkk1, kkk2, kkk3]
%for kkk4 = 1:length(LV)
    
        L1 = LV(kkk1);
        L3 = LV(kkk2); 
        %L5 = LV(kkk3); 
        %L7 = LV(kkk4); 
        xk = [L1; L3];
        
         den = [1 ;
                a ;
                (w1^2 + w3^2 + b);
                (a*w1^2 + a*w3^2 + L1*b + L3*b);
                (w1^2*w3^2 + b*w1^2 + b*w3^2);
                (a*w1^2*w3^2 + L3*b*w1^2 + L1*b*w3^2) ;
                b*w1^2*w3^2];

       
        rot = roots(den);
        rot(1) = rot(1) + 800 - 5000*i;
        rot(2) = rot(2) + 800 + 5000*i;
        plot(rot,'g.')
        
end
end

 
 


%L1 = opt_gains(1);
%L3 = opt_gains(2); 
%L5 = opt_gains(3); 
%L7 = opt_gains(4); 
       
L1 = 1;
L3 = 1;
        
 den = [1 ;
                a ;
                (w1^2 + w3^2 + b);
                (a*w1^2 + a*w3^2 + L1*b + L3*b);
                (w1^2*w3^2 + b*w1^2 + b*w3^2);
                (a*w1^2*w3^2 + L3*b*w1^2 + L1*b*w3^2) ;
                b*w1^2*w3^2];
            
        rot = roots(den);
        rot(1) = rot(1) + 800 - 5000*i;
        rot(2) = rot(2) + 800 + 5000*i;
        %plot(rot,'ro','MarkerSize',3)
        plot(rot,'bo')
        
        p11 = [real(rot(1)) imag(rot(1))];                         % First Point
        p21 = [real(rot(3)) imag(rot(3))];    
        p31 = [real(rot(5)) imag(rot(5))];    
L1 = 500;
L3 = 200;
        
 den = [1 ;
                a ;
                (w1^2 + w3^2 + b);
                (a*w1^2 + a*w3^2 + L1*b + L3*b);
                (w1^2*w3^2 + b*w1^2 + b*w3^2);
                (a*w1^2*w3^2 + L3*b*w1^2 + L1*b*w3^2) ;
                b*w1^2*w3^2];
            
        rot = roots(den);
        rot(1) = rot(1) + 800 - 5000*i;
        rot(2) = rot(2) + 800 + 5000*i;
        %plot(rot,'ro','MarkerSize',3)
        plot(rot,'ro')
        p12 = [real(rot(1)) imag(rot(1))]; 
        p22 = [real(rot(3)) imag(rot(3))];
        p32 = [real(rot(5)) imag(rot(5))];
        
fontsize = 7;
set( 0, 'defaultAxesFontSize', fontsize);
        
 %plot([0,real(rot(1))],[0,imag(rot(1))],'r')
 %plot([0,real(rot(2))],[0,imag(rot(2))],'r')
 %axis([-1200,100,-2500,2500])

       
 
 dp1 = p12-p11;                         % Difference
 dp2 = p22-p21; 
 dp3 = p32-p31; 

 quiver(p11(1),p11(2),dp1(1),dp1(2),0, 'MaxHeadSize',0.5)
 quiver(p21(1),p21(2),dp2(1),dp2(2),0, 'MaxHeadSize',0.5)
 quiver(p31(1),p31(2),dp3(1),dp3(2),0, 'MaxHeadSize',0.5)

 
 axis([-1200,100,-2000,2000])
    
 plot([0,real(rot(3))],[0,imag(rot(3))],'r')
 txt = '\alpha_i';
text(20,20,txt)

alpha = atan(-real(rot(3))/imag(rot(3)));
thvec = 0:0.01:alpha;
V1 = 150*cos(thvec + pi/2);
V2 = 150*sin(thvec+ pi/2);
plot(V1,V2,'r')

%figuresize(8.2,5.6,'cm')
print -depsc root_locus_2D
        
        
