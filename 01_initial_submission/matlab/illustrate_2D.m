close all
clear

% figure options
myFontSize = 8;
myLineWidth = 2;
%figSize = [530    55   360   480];
figSize = [530    55   360   400];

FB=input('Store as eps file (y/N)?: ','s');

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

% -----------------------------------------------------
% least squares fitting over grid

L = 0;      % lower bound of gains
U = 1000;   % upper bound of gains

h1 = 125; % step size dimension 1
h2 = 10;  % step size dimension 2
nPoints  = ((U-L)/h1);
nPoints2 = ((U-L)/h2);
AllRoots = zeros(nPoints, nPoints2, 6);

% enumerate the boxes
LV  = L:h1:U-h1; % vector of gains
LV2 = L:h2:U-h2; % vector of gains
figure(1), grid on
set(gcf,'outerposition', figSize, 'PaperPositionMode', 'auto')

% plot axes
plot([0, 0], [-2500, 2500], 'k'), hold on
plot([-1900, 300], [0, 0], 'k')

for kkk1 = 1:length(LV)
    for kkk2 = 1:length(LV2)
        
        L1 = LV(kkk1);
        L3 = LV2(kkk2);
        
        % coefficients of denominator polynomial
        den = [1 ;
            a ;
            (w1^2 + w3^2 + b);
            (a*w1^2 + a*w3^2 + L1*b + L3*b);
            (w1^2*w3^2 + b*w1^2 + b*w3^2);
            (a*w1^2*w3^2 + L3*b*w1^2 + L1*b*w3^2) ;
            b*w1^2*w3^2];
        
        rot = roots(den);
        
        % shift one pole pair closer to origin to have better view
        rot(1) = rot(1) + 800 - 5000i;
        rot(2) = rot(2) + 800 + 5000i;
        
        AllRoots( kkk1, kkk2, : ) = rot; 
    end
end

%AllRoots = complex(AllRoots);
plot(reshape(AllRoots, [6*nPoints*nPoints2, 1]), 'g.'), hold on, grid on

for ii = 1: nPoints
    for kk = 1:6
        plot(AllRoots(ii,:,kk), '.', 'color', [0+(ii-1)/(nPoints-1)*0.5, 0.5+(ii-1)/(nPoints-1)*0.5, 0])
    end
end


%%% Plot blue circles for small PR gains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
plot(rot, 'bo')

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
plot(rot, 'ro')

% plot connection to zero
plot([0,real(rot(1))], [0,imag(rot(1))], 'r')
plot([0,real(rot(3))], [0,imag(rot(3))], 'r')
axis([-1200, 100, -1500, 1500])
axis equal

%%% Plot connection between small and larger gain points %%%%%%%%%%%%%%%%%%
% p12 = [real(rot(1)) imag(rot(1))];
% p22 = [real(rot(3)) imag(rot(3))];
% p32 = [real(rot(5)) imag(rot(5))];

% dp1 = p12-p11;                         % Difference
% dp2 = p22-p21;
% dp3 = p32-p31;

% quiver(p11(1),p11(2),dp1(1),dp1(2),0, 'MaxHeadSize',0.5)
% quiver(p21(1),p21(2),dp2(1),dp2(2),0, 'MaxHeadSize',0.5)
% quiver(p31(1),p31(2),dp3(1),dp3(2),0, 'MaxHeadSize',0.5)


plot([0, real(rot(5))], [0, imag(rot(5))], 'r')
txt = '$\alpha_i$';
text(20, 100, txt, 'FontSize', myFontSize)

alpha = atan(-real(rot(1))/imag(rot(1)));
thvec = 0:0.01:alpha;
V1 = 150*cos(thvec + pi/2);
V2 = 150*sin(thvec + pi/2);
plot(V1, V2, 'r')

xlabel('Real Axis')
ylabel('Imaginary Axis')
set(gca,'FontSize', myFontSize);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)


if (FB=='y' | FB=='Y')
    matlabfrag('root_locus_2D')
    movefile('root_locus_2D.*', '../fig', 'f')
end



