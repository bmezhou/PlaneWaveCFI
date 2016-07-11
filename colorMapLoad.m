function [ map ] = colorMapLoad( ~ )
%COLORMAPLOAD Summary of this function goes here
%   Detailed explanation goes here

data = imread('E:\phd.works\2.Ultrasonix\CProgramme_607\sdk_6.0.7\porta\dat\map\Velocity\7.png');

data = data(:, 1, :);
data = squeeze(data);
% disp(data);
% 
map = double(data)/255;
% map = uint8(data);

% fp = fopen('map.dat', 'w');
% fwrite(fp, map, 'uint8');
% fclose(fp);

% a = zeros(500, 256);
% 
% % for i = 1:256
%     for j = 1:256
%         a(:, j) = j-1;
%     end
% % end
% 
% image(a);
% colormap(map);
% axis('image');

end