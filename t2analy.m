%% T2-star is obtained by fitting exponential curves on signal intensities
%%  Prepare MATLAB space
clc
clear all;
close all;

patient = '13'; % patient
field = '3.0T'; % field strength
scan = 1;   % scanning number
sec_num = 1;   % select analyze sector
path = ['I:\Iron Loading Assessment\data\' patient '\' field '\Scan ' num2str(scan) '\'];
load([path 'matfile.mat']);
temp = ones(size(Image(:,:,1)));

set(gcf,'position',[500 500 size(Image,2) size(Image,1)]);
imshow(temp,'border','tight','initialmagnification','fit');
hold on;axis off;
axis image
[c0,h0]=contour(u_outer,[0 0],'r');
[c1,h1]=contour(u_inner,[0 0],'r');
L = length(c0(1,:));
x0 = c0(1,2:L);
y0 = c0(2,2:L);
w1 = fill(x0,y0,'black');
L1 = length(c1(1,:));
x1 = c1(1,2:L1);
y1 = c1(2,2:L1);
w2 = fill(x1,y1,'w');
f=getframe(gcf);
imwrite(f.cdata,[path 'area.bmp'])
close all;

%% convert area image into binary image
clc
A = imread([path 'area.bmp']);
t = graythresh(A);
B = im2bw(A,t); %myocardium
C = B - bwareaopen(B,7000); %blood pool

%% creating segmentation sectors
clc
x0 = x_0;
y0 = y_0;
r = 50;
h = [];
i = 0;
for theta = 0:2/3*pi:2*pi
    x = x0 + r*cos(theta);
    y = y0 + r*sin(theta);
    i = i+1;
    h(i) = line([x0 x],[y0, y]);  
end
t = i;  
%  set(h(2),'color','r'); %h(1)不是第一条线的句柄，
%   h（13）为第一条直线句柄   ？？？
T = {};
j = 1;
for i=2:t-1
    T{j} = getAssMatrix(h(i),h(i+1),Image(:,:,1));
    j = j+1;
end
T{j} = getAssMatrix(h(i),h(i+1),Image(:,:,1));
close all;
%% Matrix for small sector segmentation
sector_myo = {};
sector_blood = {};
si_myo = [];
si_blood = [];
Area_myo = [];
Area_blood = [];
for i = 1:j
    sector_myo{i} = and(~B,T{i});
    sector_blood{i} = and(C,T{i});
    % image after alignment and noise remove
    [si_myo(i,:),Area_myo(i,:)] = getImageIntensity(sector_myo{i},Image(:,:,:));
    [si_blood(i,:),Area_blood(i,:)] = getImageIntensity(sector_blood{i},Image(:,:,:));
end

% save matfile\si.mat si

%% 计算出回波弛豫TE的值
te = [];
for i=1:8
    te(i) = TE_TR(i,scan);
end

%% Obtain T2* values of myocardium sectors by exponential curve fitting
clc

a_myo = [];
t2star_myo = [];
for i= 1:j
    f = @(p,te)p(1)*exp(p(2)*te);
    x0 = [100,-1];
    p_myo = lsqcurvefit(f,x0,te,si_myo(i,:));
    a_myo(i) = p_myo(1);
    t2star_myo(i) = -1/p_myo(2);
end
si_myo(:,9) = t2star_myo;
f1 = figure(1);
set(f1,'position',[33 23 1191 500]);
t = uitable;
rnames = {};
 for i=1:8
    rnames{i}=['TE' num2str(TE_TR(i,scan))];   
 end
rnames{9}='T2*';

 set(t,'RowName',rnames,'Data',si_myo');
 set(t,'Position',[100 100 1000 300]);

%% T2* mapping
clc
figure(2);
set(gcf,'Position',[32 502 989 420]);
subplot(121);
imagesc(Image(:,:,1));colormap(gray);
axis off;
hold on;
x0 = x_0;
y0 = y_0;
r = 50;
for theta = 0:2/3*pi:2*pi
    x = x0 + r*cos(theta);
    y = y0 + r*sin(theta);
    line([x0 x],[y0, y],'color','r');  
end
contour(u_outer,[0 0],'c');
contour(u_inner,[0 0],'c');

c=contour(sector_myo{sec_num},[0 0]);
L = length(c(1,:));
fill(c(1,2:L),c(2,2:L),'g'); 
subplot(122);
SI_myo = si_myo(sec_num,1:8);
t2 = t2star_myo(sec_num);
b = a_myo(sec_num);
T2curve(t2,te,SI_myo,b);

%% Obtain T2* values of blood pool sectors by exponential curve fitting
clc

a_blood = [];
t2star_blood = [];
for i= 1:j
    f = @(p,te)p(1)*exp(p(2)*te);
    x0 = [100,-1];
    p_blood = lsqcurvefit(f,x0,te,si_blood(i,:));
    a_blood(i) = p_blood(1);
    t2star_blood(i) = -1/p_blood(2);
end
si_blood(:,9) = t2star_blood;
f2 = figure(3);
set(f2,'position',[33 23 1191 500]);
t = uitable;
rnames = {};
 for i=1:8
    rnames{i}=['TE' num2str(TE_TR(i,scan))];   
 end
rnames{9}='T2*';

 set(t,'RowName',rnames,'Data',si_blood');
 set(t,'Position',[100 100 1000 300]);

%% T2* mapping
clc
figure;
set(gcf,'Position',[32 502 989 420]);
subplot(121);
imagesc(Image(:,:,1));colormap(gray);
axis off;
hold on;
x0 = x_0;
y0 = y_0;
r = 50;
for theta = 0:2/3*pi:2*pi
    x = x0 + r*cos(theta);
    y = y0 + r*sin(theta);
    line([x0 x],[y0, y],'color','r');  
end
contour(u_outer,[0 0],'c');
contour(u_inner,[0 0],'c');

c=contour(sector_blood{sec_num},[0 0]);
L = length(c(1,:));
fill(c(1,2:L),c(2,2:L),'g'); 
subplot(122);
SI_blood = si_blood(sec_num,1:8);
t2 = t2star_blood(sec_num);
b = a_blood(sec_num);
T2curve(t2,te,SI_blood,b);