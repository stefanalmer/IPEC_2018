close all
clear

% figure options
myFontSize = 8;
myLineWidth = 2;
figSize = [530    55   360   480];
FB=input('Store as eps file (y/N)?: ','s');

%yalmip('clear')
%oldopts = sdpsettings;
%solveropts = sdpsettings(oldopts,'verbose',0,'solver','sedumi');


% options = optimset;
% options = optimset(options,'Display','off');

% transfer function
L    =  85e-6;   % [H]
Ccap = 275e-6;   % [F]
R    = 0.010;    % [ohm]

% G = b/(s^2 + a*s + b)
b = 10463/(Ccap*10000*L);
a = (3259 + 10000*R)/(10000*L);

% harmonics
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

L = 0;      % lower bound of gains
U = 1000;   % upper bound of gains

% C = [-1; 0; 0; 0; 0];   % grid size of boxes
h1 = 125; % 100; % step size
h2 = 10;
nPoints = ((U-L)/h1);
nPoints2 = ((U-L)/h2);
AllRoots = zeros(nPoints, nPoints2, 6);

% enumerate the boxes
LV = L:h1:U-h1; % vector of gains
LV2 = L:h2:U-h2; % vector of gains
figure(1), grid on
set(gcf,'outerposition', figSize, 'PaperPositionMode', 'auto')

% plot axes
plot([0, 0], [-2500, 2500], 'k'), hold on
plot([-1900, 300], [0, 0], 'k')

for kkk1 = 1:length(LV)
    for kkk2 = 1:length(LV2)
        %for kkk3 = 1:length(LV)
        %     [kkk1, kkk2, kkk3]
        %for kkk4 = 1:length(LV)
        
        L1 = LV(kkk1);
        L3 = LV2(kkk2);
        %L5 = LV(kkk3);
        %L7 = LV(kkk4);
        %xk = [L1; L3];
        
        % coefficients of denominator polynomial
        den = [1 ;
            a ;
            (w1^2 + w3^2 + b);
            (a*w1^2 + a*w3^2 + L1*b + L3*b);
            (w1^2*w3^2 + b*w1^2 + b*w3^2);
            (a*w1^2*w3^2 + L3*b*w1^2 + L1*b*w3^2) ;
            b*w1^2*w3^2];
        
        
        rot = roots(den);
        rot(1) = rot(1) + 800 - 5000i;
        rot(2) = rot(2) + 800 + 5000i;
        
        AllRoots( kkk1, kkk2, : ) = rot; 
        %plot(rot,'g.'), hold on
        
    end
end
AllRoots = complex(AllRoots);
plot(reshape(AllRoots, [6*nPoints*nPoints2, 1]), 'g.'), hold on, grid on

for ii = 1: nPoints
    for kk = 1:6
        plot(complex(AllRoots(ii,:,kk)), '.', 'color', [0+(ii-1)/(nPoints-1)*0.5, 0.5+(ii-1)/(nPoints-1)*0.5, 0])
        %plot3(real(AllRoots(ii,:,kk)), imag(AllRoots(ii,:,kk)), 1:nPoints, '.', 'color', [0, 0, 0.2+ii/nPoints*0.8])
    end
end



%%% Plot blue circles for small PR gains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
rot(1) = rot(1) + 800 - 5000i;
rot(2) = rot(2) + 800 + 5000i;
%plot(rot,'ro','MarkerSize',3)
plot(rot,'bo')

p11 = [real(rot(1)) imag(rot(1))];                         % First Point
p21 = [real(rot(3)) imag(rot(3))];
p31 = [real(rot(5)) imag(rot(5))];


%%% Plot red circles for larger PR gains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
rot(1) = rot(1) + 800 - 5000i;
rot(2) = rot(2) + 800 + 5000i;
%plot(rot,'ro','MarkerSize',3)
plot(rot,'ro')
p12 = [real(rot(1)) imag(rot(1))];
p22 = [real(rot(3)) imag(rot(3))];
p32 = [real(rot(5)) imag(rot(5))];

% fontsize = 7;
% set( 0, 'defaultAxesFontSize', fontsize);

plot([0,real(rot(1))],[0,imag(rot(1))], 'r')
plot([0,real(rot(3))],[0,imag(rot(3))], 'r')
axis([-1200, 100, -1500, 1500])
axis equal

%%% Plot connection between small and larger gain points %%%%%%%%%%%%%%%%%%
dp1 = p12-p11;                         % Difference
dp2 = p22-p21;
dp3 = p32-p31;

% quiver(p11(1),p11(2),dp1(1),dp1(2),0, 'MaxHeadSize',0.5)
% quiver(p21(1),p21(2),dp2(1),dp2(2),0, 'MaxHeadSize',0.5)
% quiver(p31(1),p31(2),dp3(1),dp3(2),0, 'MaxHeadSize',0.5)




plot([0,real(rot(5))],[0,imag(rot(5))], 'r')
txt = '$\alpha_i$';
text(20,100, txt, 'FontSize', myFontSize)

alpha = atan(-real(rot(1))/imag(rot(1)));
thvec = 0:0.01:alpha;
V1 = 150*cos(thvec + pi/2);
V2 = 150*sin(thvec+ pi/2);
plot(V1,V2,'r')

xlabel('Real')
ylabel('Imaginary')
set(gca,'FontSize', myFontSize);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)
%figuresize(8.2,5.6,'cm')

if (FB=='y' | FB=='Y')
    %print -depsc root_locus_2D
    matlabfrag('root_locus_2D')
    movefile('root_locus_2D.*', '../fig', 'f')
end



