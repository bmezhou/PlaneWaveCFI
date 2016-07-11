% close all
clear all
clc

% parameters
nCh = 128;              % number of channels
reRoute = true;         % true: transducer element (correct image), false: DAQ element
chanls = ones(1, nCh);  % what channels to read (DAQ element), for each channel set to 1 
                        % if you want to read the data

% Folder path
path = 'C:\DaqData\20160518'; %Carotid  518
if (path(end) ~= '\') 
    path = [path,'\'];
end

[hdr, RF] = readDAQ(path, chanls, 1, reRoute);
bfD = zeros(size(RF,1), size(RF,2) * 2);

% bfDas = zeros(hdr(3), 257, hdr(2));
%% Plane-wave imaging

[ape, delay, apeWin] = beamFormMakeTab(RF, 73, 0);

[~, RF] = readDAQ(path, chanls, 1, reRoute);

tic;
for i = 1:1
    bfDas0 = beamFormRect(RF, 73, 0);
end
t = toc;
disp(t);

tic;
for i = 1:1
    bfDas1 = beamFormLUT(RF, 73, 0, ape, delay, apeWin);
end
t = toc;
disp(t);

tic;
for i = 1:1
    bfDas2 = beamFormLUT(RF, 73, 0, ape, delay, apeWin);
end
t = toc;
disp(t);

% for i = 1:10: hdr(2)
%     [~, RF] = readDAQ(path, chanls, i, reRoute);
%     [bfDas(:,:,i), ~]
%     bfDas0 = beamFormRect(RF, 73, 0);
%     bfDas1 = beamFormLUT(RF, 73, 0, ape, delay);


% for i = 1:10:hdr(2)
    bfDas = bfDas1(200: 1800, :);
    env = abs(hilbert(bfDas));
    logEnv = (20*log10(env/max(max(env))) + 50 ) * 255/80;

    y = 0: 1: size(bfDas, 1) - 1;
    y = y / 40e6 * 1540 / 2;

    x = 0: 1: size(bfDas, 2) - 1;
    x = x * 0.3048e-3 / 2;

    [xo, yo] = meshgrid(x, y);
    [xi, yi] = meshgrid(0:0.1e-3:max(x), 0:0.1e-3:max(y));

    logEnvIntp = interp2(xo, yo, logEnv, xi, yi);

    xn = 0:0.1e-3:max(x);
    yn = 0:0.1e-3:max(y);

    figure;
%     image(xn, yn, logEnvIntp);
    image(logEnvIntp);
    colormap(gray(256));
    axis('image');
    title(i)
    drawnow;

    
    h = imrect(gca, [10 10 100 100]);
%     h = imrect(gca, [0.005 0.005 0.0100 0.0100]);
    addNewPositionCallback(h, @(p) save('pos.mat', 'p'));
    fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(h,fcn); 

%}
%%
%{
envOri = abs(hilbert(RF));
logOri = (20*log10(envOri/max(max(envOri))) + 80 ) * 255/80;

figure;
image(logOri);
colormap(gray(256));
%%
figure;
for i = 1: 18
    plot(RF(:, i) + (i - 2) * max(max(RF)), 'r');
    hold on
end
%}