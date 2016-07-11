function [logEnvIntp] = bModeImg(ape, delay, apeWin, maxEnv, frameNum)

% parameters
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
% [hdr, RF] = readDAQ(path, chanls, 1, reRoute);
% bfD = zeros(size(RF,1), size(RF,2) * 2);
% 
% % bfDas = zeros(hdr(3), 257, hdr(2));
% %% Plane-wave imaging
% tic;
% [ape, delay, apeWin] = beamFormMakeTab(RF, 73, 0);
% disp(toc);
[~, RF] = readDAQ(path, chanls, frameNum, reRoute);

% tic;
bfDas1 = beamFormLUT(RF, 73, 0, ape, delay, apeWin);
% disp(toc);

% for i = 1:10: hdr(2)
%     [~, RF] = readDAQ(path, chanls, i, reRoute);
%     [bfDas(:,:,i), ~]
%     bfDas0 = beamFormRect(RF, 73, 0);
%     bfDas1 = beamFormLUT(RF, 73, 0, ape, delay);


% for i = 1:10:hdr(2)
bfDas = bfDas1(200: 1800, :);
% tic;
env = abs(hilbert(bfDas));
% disp(toc);
logEnv = (20*log10(env/maxEnv) + 35 ) * 255/80;

y = 0: 1: size(bfDas, 1) - 1;
y = y / 40e6 * 1540 / 2;

x = 0: 1: size(bfDas, 2) - 1;
x = x * 0.3048e-3 / 2;

[xo, yo] = meshgrid(x, y);
[xi, yi] = meshgrid(0:0.1e-3:max(x), 0:0.1e-3:max(y));

logEnvIntp = interp2(xo, yo, logEnv, xi, yi);

xn = 0:0.1e-3:max(x);
yn = 0:0.1e-3:max(y);

logEnvIntp(logEnvIntp<0) = 0;
logEnvIntp(logEnvIntp>255) = 255;

logEnvIntp = single(logEnvIntp);

% figure;
% % image(xn, yn, logEnvIntp);
% image(logEnvIntp);
% colormap(gray(256));
% axis('image');    
end