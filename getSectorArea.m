%获取划分小扇形区域
function [sector] = getSectorArea(u,a,b)
%
Temp = zeros(size(u));
s = size(u);
f = figure;
imshow(Temp);hold on;
x = [a(1) a(2) b(1)];
y = [a(3) a(4) b(2)];
fill(x,y,'w');
saveas(gcf,'area.bmp');   %保存figure上的对象
close(f);
A = imread('area.bmp');
% subplot(221);imshow(A);
B = ~A;          %取反
% subplot(222);imshow(B);
[ym xm] = find(B,1,'first');    %找到矩阵中第一个为1的位置
B = imcrop(A,[xm ym s(2)-1 s(1)-1]);   
B = im2bw(B);
% subplot(223);imshow(B);
sector = and(u,B); %获得两矩阵相交的扇形区域
% subplot(224);imshow(sector);

