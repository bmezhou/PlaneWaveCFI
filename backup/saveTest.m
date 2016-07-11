% clear
% % clc
% 
% load C:\DaqData\20160518\anaData\analyDat.mat;
% 
% 
% dataC = data;
% 
% data = real(dataC);
% 
% tic;
% saveC(dataC);
% disp(toc);
% 
% 
% fp = fopen('data1.dat', 'wb+');
% tic;
% fwrite(fp, real(dataC), 'single');
% fwrite(fp, imag(dataC), 'single');
% disp(toc);
% fclose(fp);
% 
% fp = fopen('data1.dat','rb');
% dataR = fread(fp, 644 * 125 * 200, 'single=>single');
% dataI = fread(fp, 644 * 125 * 200, 'single');
% fclose(fp);

clear;

a = zeros(1, 100, 'single');
b = ones(1, 10, 'single');

c = sseFilter(a, b);

load bfDas.mat;

load Num1.mat;
Num = Num1;

data = bfDas(:, 128);

data = data';

% Num = 1.5;

tic;
for i = 1:2000
    y0 = filter(Num, 1, data);
end
disp(toc);

tic;
for i = 1:2000
    y1 = sseFilter(data, single(Num));
end
disp(toc);

plot(y0(length(Num)-0:end));
hold on
plot(y1, 'r');

% hold on
% plot(data, 'k');
% xlim([1, 500])

% tic
% dataR = reshape(dataR, 644, 125, 200);
% dataI = reshape(dataI, 644, 125, 200);
% % disp(toc);
% 
% err1 = sum(sum(sum(dataR - real(dataC), 1), 2), 3);
% err2 = sum(sum(sum(dataI - imag(dataC), 1), 2), 3);
% 
% tic;
% save('data0.mat', 'dataC', '-v6');
% disp(toc);
% 
% tic;
% save('data1.mat', 'data');
% disp(toc);
% 
% tic;
% savefast('data2.mat', 'dataC');
% disp(toc);
% 
% 
% filename='data3.mat';
% [filepath, filebase, ext] = fileparts('data3.mat');
% if isempty(ext)
%     filename = fullfile(filepath, [filebase '.mat']);
% end
% 
% tic;
% h5create(filename, '/data', size(data), 'DataType', class(data));
% h5write(filename, '/data', data);
% disp(toc);