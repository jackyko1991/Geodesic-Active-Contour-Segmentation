%% read image
Img1 = imread('c:\matlab\drlse\1.bmp');
Img2 = imread('c:\matlab\drlse\2.bmp');
Img3 = imread('c:\matlab\drlse\3.bmp');
Img4 = imread('c:\matlab\drlse\4.bmp');
Img5 = imread('c:\matlab\drlse\5.bmp');
Img6 = imread('c:\matlab\drlse\6.bmp');
Img7 = imread('c:\matlab\drlse\7.bmp');
Img8 = imread('c:\matlab\drlse\8.bmp');
Img = {Img1,Img2,Img3,Img4,Img5,Img6,Img7,Img8};
% for i=1:8
%     Img{i}=double(Img{i}(:,:,1));
% end
% for i=1:8
%     subplot(2,4,i);
%     imagesc(Img{i},[0,255]);
%     colormap(gray);
% end

%% image registration
% load inputpoint;
% load basepoint;
% base_points = base_points{1};
% for i=1:7
%    tform = cp2tform(input_points{i},base_points,'linear conformal');
%    Img{i+1} = imtransform(Img{i+1},tform,...
%                          'FillValue',255,...
%                          'XData',[1 size(Img1,2)],'YData',[1 size(Img1,1)]);
% end
% for i=1:8
%     subplot(4,2,i);imagesc(Img{i},[0, 255]);colormap(gray);
% end

 %% image registration
  [Maxmi1,X1,Y1,Image1] = Mcurve_My(Img1,Img2,5,70);
  [Maxmi2,X2,Y2,Image2] = Mcurve_My(Img1,Img3,5,70);
  [Maxmi3,X3,Y3,Image3] = Mcurve_My(Img1,Img4,5,70);
  [Maxmi4,X4,Y4,Image4] = Mcurve_My(Img1,Img5,5,70);
  [Maxmi5,X5,Y5,Image5] = Mcurve_My(Img1,Img6,5,70);
  [Maxmi6,X6,Y6,Image6] = Mcurve_My(Img1,Img7,5,70);
  [Maxmi7,X7,Y7,Image7] = Mcurve_My(Img1,Img8,5,70);
  Image = {Img1,Image1,Image2,Image3,Image4,Image5,Image6,Image7};
  for i=1:8
      subplot(2,4,i);imagesc(Image{i});colormap(gray);
  end

  
 %% initial contour
    [nrow ncol] = size(Img1);
    initialLSF = sdf2circle(nrow,ncol,nrow/2,ncol/2-16,23);  %23,30
    initialLSF_0 = sdf2circle(nrow,ncol,nrow/2,ncol/2-16,15);%15,20
    u = initialLSF;
    u_0 = initialLSF_0;
    for i=1:8
         subplot(2,4,i);imagesc(Image{i});
         colormap(gray);hold on;
         [c0,h0] = contour(u,[0,0],'r');
         [c1,h1] = contour(u_0,[0,0],'r');title('Initial contour');
    end

   
%% processing image
sigma=1.2;    % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma);

Imgsmooth=cell(1,8);
f = cell(1,8);
g = cell(1,8);
for i=1:8
    Imgsmooth{i}= conv2(double(Image{i}),double(G),'same');  %图像高斯平滑去噪
    [Ix Iy] = gradient(Imgsmooth{i}); 
    f{i} = Ix.^2+Iy.^2;
    g{i} = 1./(1+f{i});    
end
% %% save the data of Ix,Iy and g{1}
% % [Ix Iy] = gradient(Imgsmooth{1}); 
% % % disp(Ix);
% % % disp(Iy);
% % A = g{1};
% % save Ix.txt -ascii Ix
% % save Iy.txt -ascii Iy
% % save g1.txt -ascii A
% 
% %% dispaly Imgsmooth
% % for i=1:8
% %     subplot(4,2,i);
% %     imagesc(Imgsmooth{i},[0,255]);
% %     colormap(gray); 
% % end
% 
% %% display f
% % for i=1:8
% %     subplot(4,2,i);
% %     imagesc(f{i},[0,255]);
% %     colormap(gray); 
% % end
% 
% %% display g
% % for i=1:8
% %     disp('the g is:');
% %     disp(g{i});
% %    
% % end
% 
% %% imshow g
% % for i=1:8
% %     subplot(4,2,i);
% % %     g{i} = 128*g{i};
% %     imagesc(g{i},[0,255]);
% %     colormap(gray); 
% % end
% 
% %% segement
alf=.5;  % coefficient of the pressure/balloon force
timestep=.1;  
iterNum_evo = 10; % number of iterations in level set evolution before re-initialization
iterNum_ri = 10; % number of iterations in re-initialization

for n=1:200
        u = MEVOLUTION_GAC(u, g , alf, timestep, iterNum_evo);
        u=reinit_SD(u, 1, 1, .5, iterNum_ri);
        u_0 = MEVOLUTION_GAC(u_0, g, alf, timestep, iterNum_evo);
        u_0 = reinit_SD(u_0,1,1,.5,iterNum_ri);      %保证曲线满足符号函数
        pause(0.05);
        if mod(n,20)==0
            for i=1:8
                subplot(2,4,i); imagesc(Image{i},[0,255]);colormap(gray);
                 hold on;
                [c2,h2] = contour(u,[0 0],'r');
                [c3,h3] = contour(u_0,[0 0],'r');
                iterNum=[num2str(n*iterNum_evo), ' iterations'];
                title(iterNum);
            end

        end
end
% %% save u and u_0
% save set.mat u u_0;
%% save set and image 
% save matfile\set.mat u u_0

%% load set
% load matfile\set

%% choose segment
clc
u1 = u;
u2 = u_0;
u1 = standard(u1);
u2 = standard(u2);
u_r = xor(u1,u2);
[SI1,Area1] = getImageIntensity(u_r,Img{1});
[SI2,Area2] = getImageIntensity(u_r,Img{2});
[SI3,Area3] = getImageIntensity(u_r,Img{3});
[SI4,Area4] = getImageIntensity(u_r,Img{4});
[SI5,Area5] = getImageIntensity(u_r,Img{5});
[SI6,Area6] = getImageIntensity(u_r,Img{6});
[SI7,Area7] = getImageIntensity(u_r,Img{7});
[SI8,Area8] = getImageIntensity(u_r,Img{8});
si = [SI1 SI2 SI3 SI4 SI5 SI6 SI7 SI8];
te = 2.6:2.02:16.74;
si

%% interesting area segment
imshow(Image{1});hold on;
[c h]=contour(u,[0 0],'r');
[c0 h0]=contour(u_0,[0 0],'r');
u1 = u;
segment(u1,6);
% h = get(gca,'Children');
% get(h,'Type');
% h1 = findobj(gca,'Type','line','-and','color','r');
% get(h1)

%% choose segment
u1 = u;
u2 = u_0;
imshow(u1)
u1 = standard(u1);
u2 = standard(u2);
u_r = xor(u1,u2);
% imshow(u_r);hold on;
% [c,h] = contour(u,[0 0],'r');
% [c0,h0] = contour(u_0,[0 0],'r');
% segment(u_r,6);
% dim = size(u_r);
% col =round(dim(2)/2);
% row = find(u_r(:,col),1);
%%

[row col] = find(u_r,1);
contour = bwtraceboundary(u_r,[row col],'N',8,4000);
plot(contour(:,2),contour(:,1),'r')

% segment(Image{1},6);

% for i=1:8
%     subplot(2,4,i); imagesc(Image{i},[0 255]);colormap(gray);
%     hold on;
%     [c h]=contour(u,[0 0],'r');
%     [c0 h0]=contour(u_0,[0 0],'r');
%     segment(Image{i},6);
% end

% chose intresting area

%% T2* analysis


%%



