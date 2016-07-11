clear

load bfDas.mat;

load Num1.mat;
Num = Num1;

data = bfDas(:, 128);

data = data';

fs = 40e6;
figure;
plot(data);
hold on
plot(abs(hilbert(data)), 'r')

t = 1/fs:1/fs:(length(Num))/fs;

fc = 4.0e6;

sinWavF = sin(2*pi*fc * t) .* Num;
cosWavF = cos(2*pi*fc * t) .* Num;

[dft, ang, frq] = spectrumAns(data, fs);

xq = filter(sinWavF, 1, [zeros(1, 5) data]);
xi = filter(cosWavF, 1, [zeros(1, 5) data]);

env = sqrt(xq.^2 + xi.^2) * 2;

hold on;

plot(env(length(Num) - 1:end), 'k');


[dftA, angA, frqA] = spectrumAns(1j * xq +  xi, fs);
[dftC, angC, frqC] = spectrumAns(hilbert(data), fs);

figure;
plot(frq, dft);
hold on
plot(frqA, dftA, 'r');
hold on
plot(frqC, dftC, 'k');

[dftF, angF, frqF] = spectrumAns(1j * [sinWavF, zeros(1, 200)] + [cosWavF, zeros(1, 200)], fs);
hold on
plot(frqF, dftF);

ylim([-80, 0.5]);    

data = data(1:30);

dataT = data(1:20);

tic;
for i = 1:800
    env = hilbert(dataT);
end
disp(toc)


tic;
for i = 1:800
    x1 = filter(cosWavF, 1, data);
    y1 = filter(sinWavF, 1, data);
end
disp(toc)

tt = 1/fs: 1/fs: length(data)/fs;
sinW = sin(2 * pi * fc * tt);
cosW = cos(2 * pi * fc * tt);

tic;
for i = 1:800
    dataX = sinW .* data;
    dataY = cosW .* data;
    x1 = filter(Num, 1, dataX);
    y1 = filter(Num, 1, dataY);
end
disp(toc)


tones = ones(1, 30);
tfilt = filter(Num, 1, tones);