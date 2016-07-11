nCh = 128;              % number of channels
reRoute = true;         % true: transducer element (correct image), false: DAQ element
chanls = ones(1, nCh);  % what channels to read (DAQ element), for each channel set to 1 
                        % if you want to read the data

path = 'C:\DaqData\20160518'; %Carotid
if (path(end) ~= '\') 
    path = [path,'\'];
end

tic

for i = 1:80
    [hdr, RF] = readDAQ(path, chanls, i, reRoute);
%     
%     filename = 'CH001.daq';
%     fid = fopen([path, filename],'r');
%     
%     fread(fid, 2500, 'int16=>float');
%     
%     fclose(fid);
end

disp(toc);