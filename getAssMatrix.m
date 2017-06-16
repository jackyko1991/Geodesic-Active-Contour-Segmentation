% function T = getAssMatrix(p0,h1,h2,img)
function T = getAssMatrix(h1,h2,img)    %构建辅助矩阵
% 获取矩阵保存形式的三角形
xdata_start1=get(h1,'xdata');   %获取当前直线的坐标
ydata_start1=get(h1,'ydata');
xdata_start2=get(h2,'xdata');
ydata_start2=get(h2,'ydata');
p0 = [xdata_start1(1) ydata_start1(1)];
p1 = [xdata_start1(2) ydata_start1(2)];
p2 = [xdata_start2(2) ydata_start2(2)];

T = zeros(size(img));
% 采用面积相等法
xmax = ceil(max(max(p0(1),p1(1)),p2(1)));
xmin = floor(min(min(p0(1),p1(1)),p2(1)));
ymax = ceil(max(max(p0(2),p1(2)),p2(2)));
ymin = floor(min(min(p0(2),p1(2)),p2(2)));
for x=xmin:1:xmax
    for y=ymax:-1:ymin
        p=[x y];
        t = IsInTriangleArea(p,p0,p1,p2);
        if t==1
            T(y,x)=1;
        end
    end
end
% save matfile\T.mat T