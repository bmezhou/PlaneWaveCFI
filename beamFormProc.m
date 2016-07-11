% clc; 
clc; clear;
load pos.mat;

x1 = p(1) * 1e-4;
z1 = p(2) * 1e-4 + 200/40e6/2 * 1540 + 1e-3;

eleWid = 0.3048e-3/2;

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

path = 'C:\DaqData\20160518'; %Carotid  518
if (path(end) ~= '\') 
    path = [path,'\'];
end
[hdr, RF] = readDAQ(path, chanls, 10, reRoute);

% alpha = - alpha;

m = round(p(3)*1e-4/2.0e-4);
n = round(p(4)*1e-4/2.0e-4);

alpha = 12 / 180 * pi;

% [ape, delay, apeWin] = angleFuncTab(RF, (73), (alpha), (x), (z), (m), (round(p(3)*1e-4/eleWid)), (round(p(4)*1e-4/deltaZ)));
[ape, apeWin, delay0, delay1] = angleFuncTabPix(73, alpha, x, z, m, n);

%%
load Num1.mat;

Num = Num1;

fs = 40e6;

t = 1/fs:1/fs:(length(Num))/fs;

fc = 4e6;

sinWavF = sin(2*pi*fc * t) .* Num;
cosWavF = cos(2*pi*fc * t) .* Num;

%%
% bfDas = angleFunc(RF, (73), (alpha), (x), (z), (m), (round(p(3)*1e-4/eleWid)), (round(p(4)*1e-4/deltaZ)));
[bfDasReR, bfDasImR] = angleFuncLUTPix(RF, m, n, ape, apeWin, delay0, delay1, single(sinWavF), single(cosWavF));
% bfDas = angleFunc(RF, (73), (alpha), (x), (z), (m), (round(p(3)*1e-4/eleWid)), (round(p(4)*1e-4/deltaZ)));
[mB,nB] = size(bfDasReR);
% 
dataRe = zeros(mB, 200, 'single');
dataIm = zeros(mB, 200, 'single');
% 
% ana1 = hilbert(bfDas);
% %{

env = zeros(m,n);

for i = 1: m
    for j = 1: n
        dataR = bfDasReR((i-1)*n*40 + (j-1)*40 + 1: (i-1)*n*40 + (j-1)*40 + 40);
        dataI = bfDasImR((i-1)*n*40 + (j-1)*40 + 1: (i-1)*n*40 + (j-1)*40 + 40);
        
        env(i, j) = sqrt(dataR(10)^2 + dataI(10)^2);
    end
end
% 
env = env/max(max(env));

logEnv = (20 * log10(env) + 30) / 60 * 255;
figure;
image(logEnv');
colormap(gray(256));
axis('image');
drawnow;


disp('Beamforming begins here');
tic;
% %%
% 
for it1 = 1:50
%     dataRe = zeros(mB, nB, 200, 'single');
%     
    for it2 = 1:200
        [~, RF] = readDAQ(path, chanls, 200 * (it1 - 1) + it2, reRoute);
  
%         bfDas = angleFunc(RF, (73), (alpha), (x), (z), (m), (round(p(3)*1vve-4/eleWid)), (round(p(4)*1e-4/deltaZ)));
%         bfDas = angleFuncLUT(RF, (73), (alpha), (x), (z), (m), (round(p(3)*1e-4/eleWid)), (round(p(4)*1e-4/deltaZ)), ape, delay, apeWin);
        [bfDasReR, bfDasImR] = angleFuncLUTPix(RF, m, n, ape, apeWin, delay0, delay1, single(sinWavF), single(cosWavF));
%         anaDat = hilbert(bfDas);
        dataRe(:,  it2) = bfDasReR;
        dataIm(:,  it2) = bfDasImR;
    end

%%
% dataR = real(anaDat);
% dataI = imag(anaDat);
    cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataP', num2str(it1), '.dat'', ''wb+'');'];
    eval(cmd);
    
    fwrite(fp, dataRe, 'single');
    fwrite(fp, dataIm, 'single');
    
    fclose(fp);
end
% 
% % save C:\DaqData\20160518\anaData\analyDat.mat data;
% % % save D:\analyDat.mat data;
disp(toc);
%}