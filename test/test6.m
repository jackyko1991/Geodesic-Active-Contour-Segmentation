clear
clc
x=0:pi/100:2*pi;
y=sin(x);
h=plot(x,y);  % h为plot线的句柄handle
set(gcf,'position',[80,100,400,600])
% 将图像设置为距屏幕左下角 [80，100]像素
% 图像大小设置为400*600像素
set(gcf,'color',[1,1,1]) % 背景色设置为白色
mkdir image 
% 在当前文件夹下新建image文件夹，如果已存在会warning，不影响运行

% ========================
saveas(gcf,['image\','test1.jpg'])

% ========================
f=getframe(gcf);
imwrite(f.cdata,['image\','test2.jpg'])

