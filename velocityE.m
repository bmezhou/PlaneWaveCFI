function [ velEst, velPow, b] = velocityE( frameN, rotateLabel, b )
%VELCOCITYE Summary of this function goes here
%   Detailed explanation goes here
nFram = 200;

highFrame = frameN * 25 + 200;
fileNum   = floor((highFrame-1)/200) + 1;

relatFrame = highFrame - (fileNum - 1) * 200;

autoProcLen = 200;

filtLen = 40;   % Do not change this number!
%%
load pos.mat

m = round(p(3)*1e-4/2.0e-4);
n = round(p(4)*1e-4/2.0e-4);

% disp(fileNum);

if relatFrame < autoProcLen
%     disp('here');
    if rotateLabel == 1
        cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataP', num2str(fileNum - 1), '.dat'', ''r'');'];
    elseif rotateLabel == 2
        cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataN', num2str(fileNum - 1), '.dat'', ''r'');'];  
    else
        cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataR', num2str(fileNum - 1), '.dat'', ''r'');'];  
    end
    eval(cmd);
    
    nReadFram = autoProcLen - relatFrame;
    fseek(fp, m * n * (nFram - nReadFram) * filtLen * 4, 'bof');
    dataRe0 = fread(fp, m * n * nReadFram * filtLen, 'float=>single');
    fseek(fp, m * n * nFram * filtLen * 4, 'bof');
    fseek(fp, m * n * (nFram - nReadFram) * filtLen * 4, 'cof');
    dataIm0 = fread(fp, m * n * nReadFram * filtLen, 'float=>single');
    
    fclose(fp);
    
    if rotateLabel == 1
        cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataP', num2str(fileNum), '.dat'', ''r'');'];
    elseif rotateLabel == 2
        cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataN', num2str(fileNum), '.dat'', ''r'');'];  
    else
        cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataR', num2str(fileNum), '.dat'', ''r'');'];  
    end
    eval(cmd);
    
    dataRe1 = fread(fp, m * n * relatFrame * filtLen, 'float=>single');
    fseek(fp, m * n * nFram * filtLen * 4, 'bof');
    dataIm1 = fread(fp, m * n * relatFrame * filtLen, 'float=>single');
    fclose(fp);
    
    dataRe = [dataRe0; dataRe1];
    dataIm = [dataIm0; dataIm1];
    
else
    
    if rotateLabel == 1
        cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataP', num2str(fileNum), '.dat'', ''r'');'];
    elseif rotateLabel == 2
        cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataN', num2str(fileNum), '.dat'', ''r'');']; 
    else
        cmd = ['fp = fopen(''E:\DAQData\anaData\anaDataR', num2str(fileNum), '.dat'', ''r'');'];
    end
    
    eval(cmd); 
    
    nReadFram = autoProcLen;
    fseek(fp, m * n * (relatFrame - nReadFram) * filtLen * 4, 'bof');
    dataRe = fread(fp, m * n * nReadFram * filtLen, 'float=>single');
    fseek(fp, m * n * nFram * filtLen * 4, 'bof');
    fseek(fp, m * n * (relatFrame - nReadFram) * filtLen * 4, 'cof');
    dataIm = fread(fp, m * n * nReadFram * filtLen, 'float=>single');

    fclose(fp);   
    
end

nReadFram = autoProcLen;
%%
dataRe = reshape(dataRe, m * n * filtLen, nReadFram);   % m*n by 1 column into m by n matrix, for wall-filtering.
dataIm = reshape(dataIm, m * n * filtLen, nReadFram);

% Lower th
% dataRe = dataRe(:, 1:2:end);
% dataRe = dataIm(:, 1:2:end);
%%

downSampRate = 2; 

dataRe = dataRe(:, 2:downSampRate:end);
dataIm = dataIm(:, 2:downSampRate:end);

%%
wallFiltX = filter(b, 1, dataRe(:, 1:end), [], 2);
wallFiltY = filter(b, 1, dataIm(:, 1:end), [], 2);

%%
% [veloEstX, veloEstY, jLabel] = auto1Lag(wallFiltX(:,:,1:60), wallFiltY(:,:,1:60), 44, 10, stepSize);
[veloEstX, veloEstY] = auto1LagPix(m, n, wallFiltX, wallFiltY, length(b), filtLen);
[velPow,   ~       ] = auto0LagPix(m, n, wallFiltX, wallFiltY, length(b), filtLen);

velEst  = atan2(veloEstY, veloEstX);
end

