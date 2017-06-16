load matfile\Dimage

%%initial contour
clc
tic
[nrow ncol] = size(Image{1});
initialLSF = sdf2circle(nrow,ncol,nrow/2-1,ncol/2-16,25);  %23,30
initialLSF_0 = sdf2circle(nrow,ncol,nrow/2-1,ncol/2-16,16);%15,20
u = initialLSF;
u_0 = initialLSF_0;
% imagesc(Image{1});
% colormap(gray);
% axis off;axis image
% hold on;
% contour(u,[0,0],'r');
% contour(u_0,[0,0],'r');
% title('Initial contour of image 1');
% %%     
% processing image
sigma=1.2;    % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma);
g = {};
  for i = 1:8
    Imgsmooth= conv2(double(Image{i}),double(G),'same');  %图像高斯平滑去噪
    [Ix Iy] = gradient(Imgsmooth);
    f = Ix.^2+Iy.^2;
    g{i} = 1./(1+f);
 end
% load matfile\infogaf
% segement
alf=.5;  % coefficient of the pressure/balloon force
timestep=.1;    %时间步长
iterNum_evo = 10; % number of iterations in level set evolution before re-initialization
iterNum_ri = 10; % number of iterations in re-initialization
% u = {};
% u_0 = {};
 for i = 1:8
    for n=1:100
        %         u = evolution(u, g, alf, timestep, iterNum_evo);
        u = evolution(u, g{i}, alf, timestep, iterNum_evo);
        u = reinitBSF(u, 1, 1, .5, iterNum_ri);
%         u_0 = evolution(u_0, g, alf, timestep, iterNum_evo);
        u_0 = evolution(u_0, g{i}, alf, timestep, iterNum_evo);
        u_0 = reinitBSF(u_0,1,1,.5,iterNum_ri);
         pause(0.05);
        if mod(n,10)==0
            subplot(2,4,i);
            imagesc(Image{i});colormap(gray);axis off;axis image
            hold on;
            contour(u,[0 0],'r');
            contour(u_0,[0 0],'r');
            iterNum=[num2str(n*iterNum_evo), 'image iterations '];
            title(iterNum);
        end
    end
  end
toc
 