function [SI,Area] = getImageIntensity(u,image)
%获取感兴趣部位的图像信息强度
% u为分割后的二值图像矩阵
% image 为原始图像
s = size(u);
si = zeros(1,8);
area = zeros(1,8);
for k = 1:8
for i = 1:s(2)
    for j = 1:s(1)
        if u(j,i) == 1
            si(k)  = si(k) + double(image(j,i,k));
            area(k) = area(k) + 1;
        else
        end
    end
end
end
Area = area;
SI = si./area;   %平均信号强度
            
            