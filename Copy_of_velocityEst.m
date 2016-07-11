%%
load pos.mat

m = round(p(3)*1e-4/2.0e-4);
n = round(p(4)*1e-4/2.0e-4);

nReadFram = 200;

cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataR', num2str(2), '.dat'', ''r'');']; 

eval(cmd); 

filtLen = 40;
    
dataRe = fread(fp, m * n * nReadFram * filtLen, 'float=>single');
fseek(fp, m * n * 200 * filtLen * 4, 'bof');
dataIm = fread(fp, m * n * nReadFram * filtLen, 'float=>single');

fclose(fp);   

%%
dataRe = reshape(dataRe, m * n * filtLen, nReadFram);   % m*n by 1 column into m by n matrix, for wall-filtering.
dataIm = reshape(dataIm, m * n * filtLen, nReadFram);

% reduce the sampling rate.
downSampRate = 2;

dataRe = dataRe(:, 1:downSampRate:end);
dataIm = dataIm(:, 1:downSampRate:end);

nReadFram = size(dataRe, 2);
%%
prf = 5e3/downSampRate;

N     = 46;               % Order
Fc    = 60/prf*2;        % Cutoff Frequency % 0.05
slope = 0;                % Stopband Slope
Dstop = 3.1622776602e-05; % 0.0001;          % Stopband Attenuation
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
wallFiltX = filter(b, 1, dataRe(:, 1:nReadFram), [], 2);
wallFiltY = filter(b, 1, dataIm(:, 1:nReadFram), [], 2);

%%
% [veloEstX, veloEstY, jLabel] = auto1Lag(wallFiltX(:,:,1:60), wallFiltY(:,:,1:60), 44, 10, stepSize);
[veloEstX, veloEstY] = auto1LagPix(m, n, wallFiltX(:, 1:N + 20), wallFiltY(:, 1:N + 20), N+1, filtLen);
[velPow,   ~       ] = auto0LagPix(m, n, wallFiltX(:, 1:N + 20), wallFiltY(:, 1:N + 20), N+1, filtLen);
velEst  = atan2(veloEstY, veloEstX);

% [veloL, powL] = velocityE(15, 2);

alpha = 12/180 * pi;

% vz = (veloR + veloL)/(1 + cos(alpha));
% vx = (veloR - veloL)/sin(alpha);

% velEst = vz;%sqrt(vz.^2 + vx.^2);

%% 
% velPow = powR;

norValu = 3.0;

velEst(velEst < -1 * norValu) = 0;
velEst(velEst >= norValu) = 0;

velEst = (velEst + norValu)/2/norValu * 255;

x0 = 0:0.2e-3:(size(velEst, 2) - 1) * 0.2e-3;
y0 = 0:0.2e-3:(size(velEst, 1) - 1) * 0.2e-3;

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
nCh = 128;              % number of channels
reRoute = true;         % true: transducer element (correct image), false: DAQ element
chanls = ones(1, nCh);  % what channels to read (DAQ element), for each channel set to 1 
%                         % if you want to read the data
% 
% % Folder path
path = 'C:\DaqData\20160518'; %Carotid  518
if (path(end) ~= '\') 
    path = [path,'\'];
end
% 
[hdr, RF] = readDAQ(path, chanls, 1, reRoute);
bfD = zeros(size(RF,1), size(RF,2) * 2);


% bfDas = zeros(hdr(3), 257, hdr(2));
%% Plane-wave imaging
tic;
[ape, delay, apeWin] = beamFormMakeTab(RF, 73, 0);
%%
[env] = bModeImg(ape, delay, apeWin, 15000, 600); 

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

