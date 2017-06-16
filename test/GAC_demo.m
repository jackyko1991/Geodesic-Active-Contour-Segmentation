% this Matlab program inplements geodesic active contour model
% Copyright (c) 2004--2007 by Chunming Li
% Author: Chunming Li
% email: li_chunming@hotmail.com
% Revision by Chunming Li 3/23/2006

clear all;
close all;
Img = imread('noisyImg.bmp');  % A noisy synthetic image
Img=double(Img(:,:,1));
sigma=1.2;    % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma);
Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.

alf=.4;  % coefficient of the pressure/balloon force
timestep=.1;  
iterNum_evo = 10; % number of iterations in level set evolution before re-initialization
iterNum_ri = 10; % number of iterations in re-initialization

%define initial level set function as the signed distance to a circle
[nrow, ncol]=size(Img);  
initialLSF=sdf2circle(nrow,ncol, nrow/2,ncol/2,35);
u=initialLSF;
figure(1);
imagesc(Img);colormap(gray);hold on;
[c,h] = contour(u,[0 0],'r');  
title('Initial contour');

% start level set evolution
for n=1:200
    u=EVOLUTION_GAC(u, g , alf, timestep, iterNum_evo); % iterNum_evo iterations
    u=reinit_SD(u, 1, 1, .5, iterNum_ri);  % re-initialization
    pause(0.05);
    if mod(n,20)==0
        imagesc(Img);colormap(gray);hold on;
        [c,h] = contour(u,[0 0],'r'); 
        iterNum=[num2str(n*iterNum_evo), ' iterations'];        
        title(iterNum);
        hold off;
    end
end
imagesc(Img);colormap(gray);hold on;
[c,h] = contour(u,[0 0],'r'); 
title('Final contour');
