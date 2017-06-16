%% data Access
 prefix = 'heart_data\Image0000';
 fnum = 1:8;
 fname = [prefix num2str(fnum(1))]
 info = dicominfo(fname);
 voxel_size = [info.PixelSpacing; info.SliceThickness]'
 hWaitBar = waitbar(0,'Reading DICOM file');
 Img = cell(1,8);
 for i = 1:1:length(fnum)
     fname = [prefix num2str(fnum(i))]
     Img{i} = dicomread(fname);
     waitbar((length(fnum)-i)/length(fnum));
 end
 delete(hWaitBar);
% %%  display
% for i=1:8
%     subplot(2,4,i);imagesc(Img{i});colormap(gray);
% end

%% save image
% save matfile\Dimage.mat Img
load matfile\image
%% initial contour
    [nrow ncol] = size(Img{1});
    initialLSF = sdf2circle(nrow,ncol,nrow/2,ncol/2-17,24);  %23,30
    initialLSF_0 = sdf2circle(nrow,ncol,nrow/2,ncol/2-17,16);%15,20
    u = initialLSF;
    u_0 = initialLSF_0;
    for i=1:8
         subplot(2,4,i);imagesc(Img{i});
         colormap(gray);hold on;
         [c0,h0] = contour(u,[0,0],'r');
         [c1,h1] = contour(u_0,[0,0],'r');title('Initial contour');
    end
%     save matfile\contour1.mat u u_0
%% 转化成二阶函数
c0=2;
u = InitialBSF(u,c0);
u0 = InitialBSF(u_0,c0);
% save matfile\contour2.mat u u_0
%% processing image
sigma=1.2;    % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma);

Imgsmooth=cell(1,8);
f = cell(1,8);
g = cell(1,8);

for i=1:8
    Imgsmooth{i}= conv2(double(Img{i}),double(G),'same');  %图像高斯平滑去噪
    [Ix Iy] = gradient(Imgsmooth{i}); 
    f{i} = Ix.^2+Iy.^2;
    g{i} = 1./(1+f{i}); 

end
%%
L = length(g);
[grow gcol] = size(g{1});
gmax =g{1}; 
gmin = g{1};
% gm = zeros(size(g{1}));   
S = grow * gcol;
% 取最大值/最小值/平均值
for i=2:L
   for j = 1:S
%        
% %          gm(j) = gm(j)+ 1/8*g{i}(j);    
%         gm(j) = min(gm(j),g{i}(j));    %取最小值
        gmax(j) = max(gmax(j),g{i}(j));    %取最大值
        gmin(j) = min(gmin(j),g{i}(j));
   end
end 

 %% 使用其他二阶函数   该方法行不通
% c0 = 3;
% initialContour = c0*ones(size(Img{1}));
% initialContour(100:156,35:80)=-c0;
% u = initialContour; 
% for i=1:8
%          subplot(2,4,i);imagesc(Img{i});
%          colormap(gray);hold on;
%          [c0,h0] = contour(u,[0,0],'r');
% %          [c1,h1] = contour(u_0,[0,0],'r');title('Initial contour');
% end
load matfile\contour1.mat
%% segement
clc
tic 

alf=.5;  % coefficient of the pressure/balloon force
timestep=.1;  

iterNum_evo = 10; % number of iterations in level set evolution before re-initialization
iterNum_ri = 10; % number of iterations in re-initialization

for n=1:50
        u = GAC_0(u, gmax, alf, timestep, iterNum_evo,'single-well');
        u = reinit_SD(u, 1, 1, .5, iterNum_ri);
        u_0 = GAC_0(u_0, gmin, alf, timestep, iterNum_evo,'single-well');
        u_0 = reinit_SD(u_0,1,1,.5,iterNum_ri);      %保证曲线满足符号函数
        pause(0.05);
        if mod(n,10)==0
            for i=1:8
                subplot(2,4,i); imagesc(Img{i});colormap(gray);
                 hold on;
                [c2,h2] = contour(u,[0 0],'r');
                [c3,h3] = contour(u_0,[0 0],'r');
                iterNum=[num2str(n*iterNum_evo), ' iterations'];
                title(iterNum);
            end

        end
end
toc
% %% save u and u_0

 save matfile\Gset1.mat u u_0;

%%
% %% interesting area segment
% load set.mat
% load image.mat
% imtool(Image{1})
% Image{1}(200:205,122)
% imshow(Image{1});hold on;
% [c h]=contour(u,[0 0],'r'); hold on;
% [c0 h0]=contour(u_0,[0 0],'r');
%% save set
% save matfile\Dset.mat u u_0
%% load set
 load matfile\Dset
 load matfile\Dimage
%% choose segment
clc
u1 = u;
u2 = u_0;
u1 = standard(u1);
u2 = standard(u2);
u_r = xor(u1,u2);
imshow(u_r)

%%
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

%%
clc
p = polyfit(te,si,4);
x = 0:.1:20;
y = polyval(p,x);
plot(te,si,'o',x,y);
grid on

%% T2* analysis
clc
T2 = [];
T2(1) = 2.02/log(SI1/SI2);
T2(2) = 2.02/log(SI2/SI3);
T2(3) = 2.02/log(SI3/SI4);
T2(4) = 2.02/log(SI4/SI5);
T2(5) = 2.02/log(SI5/SI6);
T2(6) = 2.02/log(SI6/SI7);
T2(7) = 2.02/log(SI7/SI8);
te = 2.02:2.02:15
y  = 55 * exp(-te./8 ) + 550 * exp(-te./100) + ...
        0 * exp( -te./1000 )    %衰减函数
% t2b = exp(-kron( te', 1./T2 ) )
% t2 =(T2(2)*detT)/(T2(2)+detT)
t2 = T2.*y/(T2+y)


%%

% subplot(221);
% imshow(u_r);hold on;
% [c,h] = contour(u,[0 0],'r');
% [c0,h0] = contour(u_0,[0 0],'r');
% segment(u_r,6);
% get(gca)

% h = get(gca,'Children')
% get(h,'Type')
% A = findobj(gca,'Type','line','-and','Type','image')

% set(A,'color','g');

% BW = edge(u_r,'sobel','r');
% subplot(222);imshow(BW);
% axe_handle = get(gca)


%%

% u_r = bwlabel(u_r);
% D = regionprops(u_r,'area','boundingbox','centroid');
% w = [D.Area]
% NR = length(w)
% V = cat(1,D.BoundingBox)
[B, L] = bwboundaries(u_r,4,'nohole');
 rgb = label2rgb(L);      %标注图像转换成rgb图像
%  gray = rgb2gray(rgb);    %grb图像转换成灰度图像
 imshow(rgb);
 


% imshow(label2rgb(L, @jet, [.5 .5 .5]))
%%

% segment(Image{1},6);

% for i=1:8
%     subplot(2,4,i); imagesc(Image{i},[0 255]);colormap(gray);
%     hold on;
%     [c h]=contour(u,[0 0],'r');
%     [c0 h0]=contour(u_0,[0 0],'r');
%     segment(Image{i},6);
% end

% chose intresting area
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
