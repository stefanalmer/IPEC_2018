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
w5 = 5*w1;
w7 = 7*w1;

% -----------------------------------------------------
% least squares fitting over grid

L = 0;      % lower bound of gains
U = 1000;   % upper bound of gains

% grid size of boxes
h1 = 25; 	% step size dimension 1
h2 = 250;   % step size dimension 2
h3 = 250;   % step size dimension 3
h4 = 250;   % step size dimension 4
nPoints  = ((U-L)/h1);
nPoints2 = ((U-L)/h2);
nPoints3 = ((U-L)/h3);
nPoints4 = ((U-L)/h4);
AllRoots = zeros(nPoints, nPoints2, nPoints3, nPoints4, 10);

% enumerate the boxes
LV =  L:h1:U-h1; % vector of gains
LV2 = L:h2:U-h2; % vector of gains
LV3 = L:h3:U-h3; % vector of gains
LV4 = L:h4:U-h4; % vector of gains
figure(1)
set(gcf,'outerposition', figSize, 'PaperPositionMode', 'auto')

% plot axes
plot([0, 0], [-8000, 8000], 'k'), hold on, grid on
plot([-3000, 300], [0, 0], 'k')


for kkk1 = 1:length(LV)
    for kkk2 = 1:length(LV2)
        for kkk3 = 1:length(LV3)
            for kkk4 = 1:length(LV4)
                
                L1 = LV(kkk1);
                L3 = LV2(kkk2);
                L5 = LV3(kkk3);
                L7 = LV4(kkk4);
                
                % coefficients of denominator polynomial
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
                AllRoots( kkk1, kkk2, kkk3, kkk4, : ) = rot;
                %plot(rot,'g.')
                
            end
        end
    end
end



%AllRoots = complex(AllRoots);
for kk = 1:10
    ii = ceil(kk/2);
    kkRoots = reshape(AllRoots(:,:,:,:, kk), [nPoints*nPoints2*nPoints3*nPoints4, 1]);
    plot(kkRoots, '.', 'color', [0+(ii-1)/(5-1)*0.5, 0.5+(ii-1)/(5-1)*0.5, 0], 'Markersize', 4)
end



%%% Plot blue circles for small PR gains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%L1 = 1;
%L3 = 1;
%L5 = 1;
%L7 = 1;

%den = [ 1;
%    a ;
%    (w1^2 + w3^2 + w5^2 + w7^2 + b);
%    (a*w1^2 + a*w3^2 + a*w5^2 + a*w7^2 + L1*b + L3*b + L5*b + L7*b) ;
%    (w1^2*w3^2 + w1^2*w5^2 + w1^2*w7^2 + b*w1^2 + w3^2*w5^2 + w3^2*w7^2 + b*w3^2 + w5^2*w7^2 + b*w5^2 + b*w7^2);
%    (a*w1^2*w3^2 + a*w1^2*w5^2 + a*w1^2*w7^2 + a*w3^2*w5^2 + a*w3^2*w7^2 + a*w5^2*w7^2 + L1*b*w3^2 + L3*b*w1^2 + L1*b*w5^2 + L5*b*w1^2 + L1*b*w7^2 + L3*b*w5^2 + L5*b*w3^2 + L7*b*w1^2 + L3*b*w7^2 + L7*b*w3^2 + L5*b*w7^2 + L7*b*w5^2);
%    (w1^2*w3^2*w5^2 + w1^2*w3^2*w7^2 + b*w1^2*w3^2 + w1^2*w5^2*w7^2 + b*w1^2*w5^2 + b*w1^2*w7^2 + w3^2*w5^2*w7^2 + b*w3^2*w5^2 + b*w3^2*w7^2 + b*w5^2*w7^2);
%    (a*w1^2*w3^2*w5^2 + a*w1^2*w3^2*w7^2 + a*w1^2*w5^2*w7^2 + a*w3^2*w5^2*w7^2 + L1*b*w3^2*w5^2 + L3*b*w1^2*w5^2 + L5*b*w1^2*w3^2 + L1*b*w3^2*w7^2 + L3*b*w1^2*w7^2 + L7*b*w1^2*w3^2 + L1*b*w5^2*w7^2 + L5*b*w1^2*w7^2 + L7*b*w1^2*w5^2 + L3*b*w5^2*w7^2 + L5*b*w3^2*w7^2 + L7*b*w3^2*w5^2);
%    (w1^2*w3^2*w5^2*w7^2 + b*w1^2*w3^2*w5^2 + b*w1^2*w3^2*w7^2 + b*w1^2*w5^2*w7^2 + b*w3^2*w5^2*w7^2) ;
%    (a*w1^2*w3^2*w5^2*w7^2 + L7*b*w1^2*w3^2*w5^2 + L5*b*w1^2*w3^2*w7^2 + L3*b*w1^2*w5^2*w7^2 + L1*b*w3^2*w5^2*w7^2);
%    b*w1^2*w3^2*w5^2*w7^2  ];

%rot = roots(den);
%plot(rot, 'bo')


%%% Plot poles for optimal PR gains %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt_gains = [  111.1867;
    313.8555;
    344.3558;
    690.0204];

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
plot(rot, 'ro')


% plot connections to zero
plot([0,real(rot(1))], [0,imag(rot(1))], 'r')
plot([0,real(rot(2))], [0,imag(rot(2))], 'r')
axis([-2500, 100, -8000, 8000])

xlabel('Real Axis')
ylabel('Imaginary Axis')
set(gca,'FontSize', myFontSize);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', myFontSize)


if (FB=='y' | FB=='Y')
    matlabfrag('root_locus')
    movefile('root_locus.*', '../fig', 'f')
end