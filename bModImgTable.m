function [ape, delay, apeWin] = bModImgTable( frameN, rotateLabel )

nCh = 128;              % number of channels
reRoute = true;         % true: transducer element (correct image), false: DAQ element
chanls = ones(1, nCh);  % what channels to read (DAQ element), for each channel set to 1 
%                         % if you want to read the data
% 
% % Folder path
path = 'C:\DaqData\20160518'; %Carotid  518
if (path(end) ~= '\') 
    path = [path,'\'];
end
% 
[hdr, RF] = readDAQ(path, chanls, 1, reRoute);
bfD = zeros(size(RF,1), size(RF,2) * 2);


% bfDas = zeros(hdr(3), 257, hdr(2));
%% Plane-wave imaging

[ape, delay, apeWin] = beamFormMakeTab(zeros(2500, 128), 73, 0);
end
