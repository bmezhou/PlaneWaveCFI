function [ logDft, dftAng, frq ] = spectrumAns( x, fs )
%SPECTRUMPROC Summary of this function goes here
%   Detailed explanation goes here

n = length(x);
kf =           0: 1: round(n/2)-1;
kb = -floor(n/2): 1: -1;
k  = [kf kb];
frq = fftshift( k * fs/n) / 1e6;


dftX   = fftshift( abs(fft(x)) );
dftAng  = fftshift(unwrap(angle(fft(x))));
logDft = 20*log10(dftX/max(dftX));

end

