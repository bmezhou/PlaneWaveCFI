load pos.mat

m = round(p(3)*1e-4/2.0e-4);
n = round(p(4)*1e-4/2.0e-4);


load veloR.mat;
load veloL.mat;

load veloZ.mat;

load env.mat;

map = colorMapLoad();

pCorner = [p(1), p(2)];

pCorner = pCorner * 1e-4;

% pCorner(2) = pCorner(2) - 0.0012;
pCorner = round(pCorner /0.1e-3);

for frameCFM = 1:153
    velEstZ = ...
        (veloR(:,:, frameCFM) + veloL(:,:, frameCFM)) / (1+cos(12/180*pi));
        
    velEstX = ...
        (veloR(:,:, frameCFM) - veloL(:,:, frameCFM)) / sin(12/180*pi);
    
    velEstMag = sqrt(velEstZ.^2 + velEstX.^2);
    velEstAng = atan2(velEstZ, velEstX);
    
    template = ones(size(velEstAng));
    template(velEstAng>pi/2) = -1;
    template(velEstAng<-pi/2) = -1;
    
    velEstVisu = template .* velEstMag;
    
    
    velEstRe = veloZ(:,:, frameCFM);
    
%     velEst = velEstZ;
    
    velEstIntpZ = veloNor(velEstZ, pi);
    velEstIntpX = veloNor(velEstX, pi);
    
    velEstIntpR = veloNor(velEstRe, pi);
    velEstIntpV = veloNor(velEstVisu, 5*pi);
    
%     imageFusX = imageFusion(env(:,:, frameCFM)/ 255, round(velEstIntpX), map, pCorner);
%     imageFusZ = imageFusion(env(:,:, frameCFM)/ 255, round(velEstIntpZ), map, pCorner);
%     imageFusR = imageFusion(env(:,:, frameCFM)/ 255, round(velEstIntpR), map, pCorner);
    imageFusV = imageFusion(env(:,:, frameCFM)/ 255, round(velEstIntpV), map, pCorner);

% %     figure;
% subplot(2,2,1)
%     image(imageFusX);
%     axis('image');
%     title(frameCFM);
%     
% subplot(2,2,2)
%     image(imageFusZ);
%     axis('image');
%     title(frameCFM);
%     
% subplot(2,2,3)
%     image(imageFusR);
%     axis('image');
%     title(frameCFM); 
%     
% subplot(2,2,4)
    image(imageFusV);
    axis('image');
    title(frameCFM);    
    
    pause(0.01);
    drawnow;
    
end