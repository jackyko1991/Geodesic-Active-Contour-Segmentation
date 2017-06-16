%% read dicom images' imformation 
 clc
 prefix = 'heart_data\Image0000';
 fnum = 1:8;
 info = cell(1,8);
 for i=1:length(fnum)
     fname = [prefix num2str(fnum(i))];
     info{i} = dicominfo(fname); 
 end
%save matfile/info.mat
%  voxel_size = [info.PixelSpacing; info.SliceThickness]';
 hWaitBar = waitbar(0,'Reading DICOM file');
 Img = cell(1,8);
 for i = 1:1:length(fnum)
     fname = [prefix num2str(fnum(i))];
     Img{i} = dicomread(fname);
     waitbar((length(fnum)-i)/length(fnum));
 end
 delete(hWaitBar);

 %%  display
% for i=1:8
% %     Img{i} = double(Img{i});
%     subplot(2,4,i);imagesc(Img{i});colormap(gray);axis off;
% end
clc
info{8}
%% save image
save matfile\image.mat Img

%% image registration
tic
[Max1,X1,Y1,Image1] = Mcurve_My(Img{1},Img{2},5,30);
[Max2,X2,Y2,Image2] = Mcurve_My(Img{1},Img{3},5,30);
[Max3,X3,Y3,Image3] = Mcurve_My(Img{1},Img{4},5,30);
[Max4,X4,Y4,Image4] = Mcurve_My(Img{1},Img{5},5,30);
[Max5,X5,Y5,Image5] = Mcurve_My(Img{1},Img{6},5,30);
[Max6,X6,Y6,Image6] = Mcurve_My(Img{1},Img{7},5,30);
[Max7,X7,Y7,Image7] = Mcurve_My(Img{1},Img{8},5,30);
toc
%% display registrated image
Image = {Img{1},Image1,Image2,Image3,Image4,Image5,Image6,Image7};
for i= 1:8
    Image{i} = double(Image{i});
    subplot(2,4,i);imagesc(Image{i});colormap(gray);axis off;
end

%%
%   save matfile\Dimage.mat Image
load matfile\Dimage

%% initial contour
[nrow ncol] = size(Image{1});
    initialLSF = sdf2circle(nrow,ncol,nrow/2-3,ncol/2-15,24);  %23,30
    initialLSF_0 = sdf2circle(nrow,ncol,nrow/2-1,ncol/2-15,17);%15,20
    u = initialLSF;
    u_0 = initialLSF_0;
    
    for i=1:8
         subplot(2,4,i);imagesc(Image{i});axis off;
         colormap(gray);hold on;
         contour(u,[0,0],'r');
         contour(u_0,[0,0],'r');title('Initial contour');
    end
% processing image
%%
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
%     subplot(2,4,i);imagesc(g{i});colormap(gray);axis off;
end
% for i=1:8
%     subplot(2,4,i);imshow(g{i})
% end
% imtool(f{i})
save matfile\Imgsmooth.mat Imgsmooth
%%
clc
[grow gcol] = size(g{1});
gmax = g{1};
gmin = g{1};
fm = f{1};
S = grow * gcol;
% 取最大值/最小值/平均值
% T = S/2-15*256;
% T = S/2;
for i=2:8
    for j = 1:S
%         if j < T
            % %          gm(j) = gm(j)+ 1/8*g{i}(j);
%             gm(i) = max(gm(j),g{i}(j));
%             gm(j) = min(gm(j),g{i}(j));    %取最大值
%           gm(j) = gm(j)+g{i}(j);
%         else
            gmin(j) = min(gmin(j),g{i}(j));    %取最小值
            gmax(j) = max(gmax(j),g{i}(j));    %取最小值
%         end
    end
end
% gm = 1./(1+fm);
% gm = 1/2*gm;
% % % figure
% imshow(gm);
% hold on;
% contour(u,[0 0],'r');
% contour(u_0,[0 0],'r')

%% save g and f
% save matfile\infogaf.mat g f 

% load matfile\infogaf
% segement
tic
alf=0.5;  % coefficient of the pressure/balloon force 
timestep=.1;  
iterNum_evo = 10; % number of iterations in level set evolution before re-initialization
iterNum_ri = 10; % number of iterations in re-initialization

for n=1:100
         u = VEVOLUTION_GAC(u, gmax, alf, timestep, iterNum_evo);
         u = reinit_SD(u, 1, 1, .5, iterNum_ri);
         u_0 = VEVOLUTION_GAC(u_0, gmin, alf, timestep, iterNum_evo);
         u_0 = reinit_SD(u_0,1,1,.5,iterNum_ri);      %保证曲线满足符号函数
        pause(0.05);
        if mod(n,10)==0
            for i=1:8
                subplot(2,4,i); imagesc(Image{i});colormap(gray);axis off;
                 hold on;
                [c2,h2] = contour(u,[0 0],'r');
                [c3,h3] = contour(u_0,[0 0],'r');
                iterNum=[num2str(n*iterNum_evo), ' iterations'];
                title(iterNum);
            end

        end
end
toc
%% test
figure
imshow(gmax);
hold on;
contour(u,[0 0],'r');
contour(u_0,[0 0],'r')
save matfile\Dset1 u u_0 
 
x0 = 55;
y0 = 128;
p0 = [x0,y0];
j = 1;
T = {};
L = length(h);
for i=2:L-1
    T{j} = getAssMatrix(p0,h(i),h(i+1),Image{1});
    j = j+1;
end
T{j} = getAssMatrix(p0,h(2),h(i+1),Image{1});
% T = getAssMatrix(p0,h1,h2,Image{1});
u1 = u;
u2 = u_0;
u1 = standard(u1);
u2 = standard(u2);
u_r = xor(u1,u2);
sector = {};
si = [];
Area = [];
for i = 1:j
    sector{i} = and(u_r,T{i});
    [si(i,:),Area(i,:)] = getImageIntensity(sector{i},Image);
end
% sector = and(u_r,T);
% [si,Area] = getImageIntensity(sector,Image);
te = zeros(1,8);
for i=1:8
    te(i) = info{i}.EchoTime;
end

for i=1:8
    te(i) = info{i}.EchoTime;
end
t2 = [];
for k = 1:j
    for i=1:7
        t2(k,i)= (te(i+1)-te(i))/log(si(k,i)/si(k,i+1));
    end
end
% t2 = [];
% for i=1:7
%     t2(i)= (te(i+1)-te(i))/log(si(i)/si(i+1));
% end
TE = te(1):2.02:te(7);
decay = 200 * exp( -TE/20 ) + 800 * exp( -TE/80 ) + 2.4; 
y = decay;
t2star = [];
for i = 1:j
    t2star(i) = t2(i,:).*y/(t2(i,:)+y);
end
si(:,9) = t2star;
rnames = cell(1,9);
for i=1:8
    rnames{i}=['TE' num2str(info{i}.EchoTime)];   
   % rnames{i}=['TE' num2str(i)]
end
rnames{9}='T2*';
% p=subplot(2,2,1:2);
% position= get(p,'Position');
t= uitable(handles.T2Table);
% set(t,'Position',[20 200 300 200])
set(t,'RowName',rnames,'Data',si','ForegroundColor','black');
% axes(handles.axes2); 
% % imshow(T)
% SI = si(1:8);
% plot(te,SI)
% p = polyfit(te,SI,3);
% x = 0:.1:20;
% y = polyval(p,x);
% plot(te,SI,'*',x,y);
% grid on;color 'g';
% xlabel('TE','color','r');
% ylabel('SI','color','r');
% % axis([xmin,xmax,ymin,ymax]);
% axis([0,16,0,200])

% --- Executes on button press in AddSector.
