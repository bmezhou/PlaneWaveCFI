function [ velEstIntp ] = veloNor( velEst, norValu )
%VELONOR Summary of this function goes here
%   Detailed explanation goes here


velEst(velEst < -1 * norValu) = 0;
velEst(velEst >= norValu) = 0;
    
velEst = (velEst + norValu)/2/norValu * 255;

x0 = 0:0.20e-3:(size(velEst, 2) - 1) * 0.20e-3;
y0 = 0:0.20e-3:(size(velEst, 1) - 1) * 0.20e-3;

xi = 0:0.10e-3: max(x0);
yi = 0:0.10e-3: max(y0);

[x0, y0] = meshgrid(x0, y0);
[xi, yi] = meshgrid(xi, yi);

velEstIntp = interp2(x0, y0, velEst, xi, yi);

end

