load pos.mat

m = round(p(3)*1e-4/2.0e-4);
n = round(p(4)*1e-4/2.0e-4);

veloR = zeros(n, m, 220, 'single');
veloL = zeros(n, m, 220, 'single');
veloZ = zeros(n, m, 220, 'single');

[ape, delay, apeWin] = beamFormMakeTab(zeros(2500, 128), 73, 0);


env = zeros(309, 388, 220, 'single');
%%
downSampRate = 2;

prf = 5e3 / downSampRate;

N     = 80;               % Order
Fc    = 70/prf*2;        % Cutoff Frequency % 0.05
slope = 0;                % Stopband Slope
Dstop = 1e-05; % 0.0001;          % Stopband Attenuation
Dpass = 0.057501127785;   % Passband Ripple

% Calculate the coefficients using the FIRCEQRIP function.
b  = firceqrip(N, Fc, [Dstop, Dpass], 'slope', slope, 'high');
% b = single(b);
% disp('design');
% minium phase.
% b = firceqrip(N, Fc, [Dstop, Dpass], 'slope', slope, 'min', 'high');
% disp('design1');
% load Num.mat;

b = single(b);

%%
tic;
for frameCFM = 1:153
    [veloR1, ~   ] = velocityE(frameCFM, 1, b);
    [veloL1, ~   ] = velocityE(frameCFM, 2, b);
    
    [veloZ1, powR] = velocityE(frameCFM, 3, b);
    
    veloR(:, :, frameCFM) = veloR1;
    veloL(:, :, frameCFM) = veloL1;
    veloZ(:, :, frameCFM) = veloZ1;
    
    
    env(:,:, frameCFM) = ...
        bModeImg(ape, delay, apeWin, 15000, frameCFM * 25 + 125 );
end

disp(toc);

save veloR.mat veloR;
save veloL.mat veloL;
save veloZ.mat veloZ;

save env.mat   env;

%% 

velEst = veloR1;
velPow = powR;

norValu = 3.0;

velEst(velEst < -1 * norValu) = 0;
velEst(velEst >= norValu) = 0;

velEst = (velEst + norValu)/2/norValu * 255;

x0 = 0:0.20e-3:(size(velEst, 2) - 1) * 0.20e-3;
y0 = 0:0.20e-3:(size(velEst, 1) - 1) * 0.20e-3;

xi = 0:0.10e-3: max(x0);
yi = 0:0.10e-3: max(y0);

[x0, y0] = meshgrid(x0, y0);
[xi, yi] = meshgrid(xi, yi);

velEstIntp = interp2(x0, y0, velEst, xi, yi);

%%
velPow = velPow / max(max(velPow));
velPow = (20 * log10(velPow) + 60)/60 * 255;

% velPow = velPow * 255;

figure;
image(velPow);
colormap(gray(256));
axis('image');


% figure;
% image(velEstIntp);
% colormap(jet(256));
% axis('image');

% surf(double(velEst));
% view([0, 90]);
%%
%%

[env] = bModeImg(ape, delay, apeWin, 15000, frameCFM * 25 + 75 ); 

map = colorMapLoad();

load pos.mat

pCorner = [p(1), p(2)];

pCorner = pCorner * 1e-4;

% pCorner(2) = pCorner(2) - 0.0012;
pCorner = round(pCorner /0.1e-3);

imageFus = imageFusion(env/ 255, round(velEstIntp), map, pCorner);

figure;
image(imageFus);
axis('image');

