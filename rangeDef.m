% mex angleFuncLUT.c COMPFLAGS="/openmp /fp:fast $COMPFLAGS"
% mex angleFuncTab.c COMPFLAGS="/openmp  $COMPFLAGS"

% mex angleFunc.c    COMPFLAGS="/openmp /fp:fast $COMPFLAGS"

% mex angleFuncLUTPix.c COMPFLAGS="/openmp /fp:fast $COMPFLAGS"
% mex angleFuncTabPix.c COMPFLAGS="/openmp  $COMPFLAGS"
clc; clear;
load pos.mat;

x1 = p(1) * 1e-4;
z1 = p(2) * 1e-4 + 200/40e6/2 * 1540 + 1e-3;

eleWid = 0.3048e-3/2;
alpha = 12 / 180 * pi;
deltaZ = 0.10e-3;
%%
n = 1;

while (n*deltaZ < z1)
    n = n + 1;
end
%%
m = 1;

while (m*eleWid < x1)
    m = m + 1;
end
%%%%%%%%%%%%%%%%
x = m * eleWid;
z = n * deltaZ;

%%
% parameters
nCh = 128;              % number of channels
reRoute = true;         % true: transducer element (correct image), false: DAQ element
chanls = ones(1, nCh);  % what channels to read (DAQ element), for each channel set to 1 
                        % if you want to read the data

path = 'C:\DaqData\20160518'; %Carotid
if (path(end) ~= '\') 
    path = [path,'\'];
end

[hdr, RF] = readDAQ(path, chanls, 10, reRoute);

%%
m = round(p(3)*1e-4/2.0e-4);
n = round(p(4)*1e-4/2.0e-4);

[ape, apeWin, delay0, delay1] = angleFuncTabPix(73, -alpha, x, z, m, n);

load Num1.mat;

Num = Num1;

fs = 40e6;

t = 1/fs:1/fs:(length(Num))/fs;

fc = 4e6;

sinWavF = sin(2*pi*fc * t) .* Num;
cosWavF = cos(2*pi*fc * t) .* Num;

tic;
for i = 1:10
    disp(i);
    [bfDasRL, bfDasIL] = angleFuncLUTPix(RF, m, n, ape, apeWin, delay0, delay1, single(sinWavF), single(cosWavF));
end
disp(toc);
% %{
% [ape, apeWin, delay0, delay1] = angleFuncTabPix(73, -alpha, x, z, m, n);

% tic;
% for i = 1:10
%     [bfDasRR, bfDasIR] = angleFuncLUTPix1(RF, m, n, ape, apeWin, delay0, delay1, single(sinWavF), single(cosWavF));
% end
% disp(toc);

env = zeros(m,n);

filtLen = 40;

tic;
for i = 1: m
    for j = 1: n
        dataR = bfDasRL((i-1)*n*filtLen + (j-1)*filtLen + 1: (i-1)*n*filtLen + (j-1)*filtLen + filtLen);
        dataI = bfDasIL((i-1)*n*filtLen + (j-1)*filtLen + 1: (i-1)*n*filtLen + (j-1)*filtLen + filtLen);
        
        env(i, j) = sqrt(dataR(30)^2 + dataI(30)^2);
    end
end
%%
env = env/max(max(env));

logEnv = (20 * log10(env) + 30) / 60 * 255;

logEnv = logEnv';

disp(toc);
figure;
image(logEnv);
colormap(gray(256));

[ape, delay, apeWin] = beamFormMakeTab(RF, 73, 0);

axis('image');
tic;
logEnvIntp = bModeImg(ape, delay, apeWin, 4000, 10);
disp('here');
disp(toc);

map = colorMapLoad();

pCorner = [p(1), p(2)];

pCorner = pCorner * 1e-4;

%% Matched here!!!
% pCorner(2) = pCorner(2) - 1e-3; %%
%%
tic;
pCorner = round(pCorner /0.1e-3);

logEnv(logEnv <   0) = 0;
logEnv(logEnv > 255) = 255;

x0 = 0:0.2e-3:(size(logEnv, 2) - 1) * 0.2e-3;
y0 = 0:0.2e-3:(size(logEnv, 1) - 1) * 0.2e-3;

xi = 0:0.10e-3: max(x0);
yi = 0:0.10e-3: max(y0);

[x0, y0] = meshgrid(x0, y0);
[xi, yi] = meshgrid(xi, yi);

logEnv = interp2(x0, y0, logEnv, xi, yi);

disp(toc);
logEnv = round(logEnv);

tic;

imageSup = imageFusion(logEnvIntp/255, logEnv, gray(256), pCorner);

disp(toc);

figure;
image(imageSup);
axis('image');

%%
% for i = 1: m
%     for j = 1: n
%         dataR = bfDasRR((i-1)*n*30 + (j-1)*30 + 1: (i-1)*n*30 + (j-1)*30 + 30);
%         
%         [xi, xq] = demod(double(dataR), 4e6, 40e6, 'qam');
%         
%         envT = sqrt(xi.^2 + xq.^2);
% %         envT = demod(hilbert(dataR));
%         env(i, j) = envT(15);
%     end
% end
% % 
% env = env/max(max(env));
% 
% logEnv = (20 * log10(env) + 60) / 60 * 255;
% figure;
% image(logEnv');
% colormap(gray(256));
% axis('image');

% figure;
% image(logEnv);

% colormap(gray(256));

% axis('image');

% figure;
% 
% for i = 100: 100
%     
%      plot(bfDas1(:, i)/max(bfDas1(:, i))); 
%      hold on; 
%      plot(bfDas0(:, i)/max(bfDas0(:, i)), 'r');
%      
%      drawnow;
% %      hold on
% %      plot(bfDas1(:, i) - bfDas0(:, i), 'k');
%      pause(0.1)
% %      clf;
% end
%}