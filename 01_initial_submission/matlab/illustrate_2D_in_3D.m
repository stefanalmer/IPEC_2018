% plot the damping of the poles as function of the gains

close all
clear

% figure options
myFontSize = 8;
myLineWidth = 2;
figSize = [530    55   360   360];
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

h1 = 4; % step size dimension 1
h2 = 4;  % step size dimension 2
nPoints  = ((U-L)/h1);
nPoints2 = ((U-L)/h2);
dMat = zeros(nPoints, nPoints2, 3);

% enumerate the boxes
LV  = L:h1:U-h1; % vector of gains
LV2 = L:h2:U-h2; % vector of gains



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
        
        % calculate damping of each pole pair
        if any(real(rot) > 0) || any(imag(rot([1,3,5])) < 0)
            dMat( kkk1, kkk2, 1:3 ) = NaN;
        else
        for ii = 1:3
            %alpha = atan(-real(rot(2*ii-1))/imag(rot(2*ii-1)));
            dMat( kkk1, kkk2, ii ) = atan(-real(rot(2*ii-1))/imag(rot(2*ii-1))); 
        end
        end
    end
end


for kk = 1:3
    figure
    surf(LV, LV2, dMat(:,:,kk)*180/pi); 
    shading flat
end

figure
set(gcf, 'outerposition', figSize, 'PaperPositionMode', 'auto')
dMat2 = min(dMat, [], 3);
surf(LV, LV2, dMat2*180/pi);
shading flat
view([-75, 40])


xlabel('$\lambda_1$')
ylabel('$\lambda_3$')
zlabel('damping [deg]')
set(gca,'FontSize', myFontSize);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)


if (FB=='y' | FB=='Y')
    matlabfrag('root_locus_2D_in_3D', 'dpi', 300)
    movefile('root_locus_2D_in_3D.*', '../fig', 'f')
end





