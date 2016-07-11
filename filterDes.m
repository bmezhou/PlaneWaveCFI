function Hd = filterDes
%FILTERDES Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 8.1 and the DSP System Toolbox 8.4.
% Generated on: 17-Jun-2016 20:54:55

% FIR constrained equiripple Lowpass filter designed using the FIRCEQRIP
% function.

% All frequency values are normalized to 1.

N     = 7;               % Order
Fc    = 0.25;            % Cutoff Frequency
slope = 0;               % Stopband Slope
Dstop = 0.00177827941;   % Stopband Attenuation
Dpass = 0.057501127785;  % Passband Ripple

% Calculate the coefficients using the FIRCEQRIP function.
b  = firceqrip(N, Fc, [Dpass, Dstop], 'slope', slope, 'min');
Hd = dfilt.dffir(b);

% [EOF]