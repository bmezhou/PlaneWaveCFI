% clc; 
clear;

load pos.mat

x1 = p(1) * 1e-4;
z1 = p(2) * 1e-4 + 200/40e6/2 * 1540;

eleWid = 0.3048e-3/2;

alpha = 12 / 180 * pi;

deltaZ = eleWid/2/tan(abs(alpha))/20;


%%
n = 1;

while (n*deltaZ < z1)
    n = n + 1;    
end

%%
m = 1;

while (m*eleWid + n*deltaZ * tan(abs(alpha)) < x1)
    m = m + 1;
end

x = m * eleWid + n * deltaZ * tan(alpha);
z = n * deltaZ;

%%
m1 = 1;

while (abs(m*eleWid + 2 * n * deltaZ * tan(alpha) - m1 * eleWid) > eleWid /2 )
    m1 = m1 + 1;
end

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

m = round(p(3)*1e-4/1.5e-4);
n = round(p(4)*1e-4/1.5e-4);


% [ape, delay, apeWin] = angleFuncTab(RF, (73), (alpha), (x), (z), (m), (round(p(3)*1e-4/eleWid)), (round(p(4)*1e-4/deltaZ)));
[ape, apeWin, delay0, delay1] = angleFuncTabPix(73, -alpha, x, z, m, n);

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
dataRe = zeros(mB, nB, 200, 'single');
dataIm = zeros(mB, nB, 200, 'single');
% 
% ana1 = hilbert(bfDas);
% %{

env = zeros(m,n);

for i = 1: m
    for j = 1: n
        dataR = bfDasReR((i-1)*n*20 + (j-1)*20 + 1: (i-1)*n*20 + (j-1)*20 + 20);
        dataI = bfDasImR((i-1)*n*20 + (j-1)*20 + 1: (i-1)*n*20 + (j-1)*20 + 20);
        
        env(i, j) = sqrt(dataR(10)^2 + dataI(10)^2);
    end
end
% 
env = env/max(max(env));

logEnv = (20 * log10(env) + 60) / 60 * 255;
figure;
image(logEnv');
colormap(gray(256));
axis('image');

disp('hello');
tic;
% %%
% 
for it1 = 1:35
%     dataRe = zeros(mB, nB, 200, 'single');
    cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataN', num2str(it1), '.dat'', ''r'');'];
    eval(cmd);
    
    dataRe = fread(fp, m * n * 20 * 200, 'single=>single');
    dataIm = fread(fp, m * n * 20 * 200, 'single=>single');
    
    fclose(fp);
    
    for it2 = 1:50:200

%         [~, RF] = readDAQ(path, chanls, it1 , reRoute);
  
%         bfDas = angleFunc(RF, (73), (alpha), (x), (z), (m), (round(p(3)*1e-4/eleWid)), (round(p(4)*1e-4/deltaZ)));
%         bfDas = angleFuncLUT(RF, (73), (alpha), (x), (z), (m), (round(p(3)*1e-4/eleWid)), (round(p(4)*1e-4/deltaZ)), ape, delay, apeWin);
%         [bfDasReR, bfDasImR] = angleFuncLUTPix(RF, m, n, ape, apeWin, delay0, delay1, single(sinWavF), single(cosWavF));
%         anaDat = hilbert(bfDas);

        bfDasReR = dataRe((it2-1) * m * n * 20 + 1: it2 * m * n * 20);
        bfDasImR = dataIm((it2-1) * m * n * 20 + 1: it2 * m * n * 20);

        env = zeros(m,n);
       
        for i = 1: m
            for j = 1: n
                dataR = bfDasReR((i-1)*n*20 + (j-1)*20 + 1: (i-1)*n*20 + (j-1)*20 + 20);
                dataI = bfDasImR((i-1)*n*20 + (j-1)*20 + 1: (i-1)*n*20 + (j-1)*20 + 20);
        
                env(i, j) = sqrt(dataR(10)^2 + dataI(10)^2);        
            end
        end
% 
        env = env/max(max(env));

        logEnv = (20 * log10(env) + 60) / 60 * 255;
% figure;
        image(logEnv');
        colormap(gray(256));
        axis('image');
        title([num2str((it1 - 1)*200 + it2),' ', num2str(it1)]);
        drawnow;
        
        pause(0.05)
        
    end
%%
% dataR = real(anaDat);
% dataI = imag(anaDat);


end
% 
% % save C:\DaqData\20160518\anaData\analyDat.mat data;
% % % save D:\analyDat.mat data;
disp(toc);
%}