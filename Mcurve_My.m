function [Maxm,X,Y,Img] = Mcurve_My(rimage1,rimage2,step,rfield)
%   MCURVE_MY    绘制图像配准测度曲线和配准图像
%   参数说明：    image1，image2为参考图像和模板图像
%                step为计算测度函数时，在x和y方向平移时的步长
%                rfield配准区域选择
%                method选择测度函数，默认为最大互信息配准MI
%   输出参数：    Maxm最大测度函数值
%                X，Y配准时x、y方向的平移长度
%
%   作者：
%   版本：V1.0
%   修改：V1.1   世界坐标系，图像中心为原点    2009.07.23
%   修改：V1.2   配准区域参数选择
% 灰度级
if isa(rimage1,'uint16')
    image1 = double(rimage1)/65535*255;
else
    image1 = double(rimage1);
end
if isa(rimage2,'uint16')
    image2 = double(rimage2)/65535*255;
else
    image2 = double(rimage2);
end
% 参考图像及大小
image1 = round(image1);
image2 = round(image2);
[mr,nr] = size(image1);
[mf,nf] = size(image2);
% 2009.07.23世界坐标系原点——图像中心
wzx1 = (mr-1)/2;
wzy1 = (nr-1)/2;
wzx2 = (mf-1)/2;
wzy2 = (nf-1)/2;

% if method == 'MI';
%     % 计算互信息方式
%     Mi_method = questdlg('使用哪种计算MI方法','方式选择','互信息','归一化互信息','互信息');
% end
k = 1;
for x = wzx1-rfield+wzx2:step:wzx1+rfield+wzx2 % 浮动图像相对参考图像平移+-rfield像素
    l = 1;
    if x <= mr
        for y = wzy1-rfield+wzy2:step:wzy1+rfield+wzy2
            if y <= nr
                Im_F = image2(mf-x+1:mf,nf-y+1:nf);
                Im_R = image1(1:x,1:y);
            else
                Im_F = image2(mf-x+1:mf,1:nf+nr+1-y);
                Im_R = image1(1:x,y-nf+1:nr);
            end
%             if method == 'MI'
                MeV(k,l) = mi(Im_R,Im_F);
                l = l+1;
%             end
        end
    else
        for y = wzy1-rfield+wzy2:step:wzy1+rfield+wzy2
            if y <= nr
                Im_F = image2(1:mf+mr+1-x,nf-y+1:nf);
                Im_R = image1(x-mf+1:mr,1:y);
            else
                Im_F = image2(1:mf+mr+1-x,1:nf+nr+1-y);
                Im_R = image1(x-mf+1:mr,y-nf+1:nr);
            end
%             if method == 'MI'
                MeV(k,l) = mi(Im_R,Im_F);
                l = l+1;
%             end
        end
    end
    k = k+1;
end
% 输出
% x = wzx1-rfield+wzx2:step:wzx1+rfield+wzx2;
% y = wzy1-rfield+wzy2:step:wzy1+rfield+wzy2;
% x = x-(wzx1+wzx2);y = y-(wzy1+wzy2);
% [x,y] = meshgrid(x,y);
% mesh(x,y,MeV);
% title('测度曲线');

[Maxm,ind] = max(MeV(:));
[X,Y] = ind2sub(size(MeV),ind);
X = (X-1)*step;Y = (Y-1)*step; % 平移像素个数
% 模板图像分别沿x，y轴平移(X-1)*step，(Y-1)*step像素点后与参考图像配准
% 此时，浮动图像中心位于参考图像坐标系(wzx1-rfield+(X-1)*step,wzy1-rfield+(Y-1)*step)处
se = translate(strel(1),[X-rfield Y-rfield]);

Img = imdilate(rimage2,se);
% subplot(121);imagesc(image2);colormap(gray);
% subplot(122);imagesc(Img);colormap(gray);
% H.Position = [132 258 260 260];
% figure(H)
% % subplot(1,3,1);
% imagesc(image2)
% title('原始浮动图像');
% colormap(gray)
% % subplot(1,3,2);
% H.Position = [402 258 260 260];
% figure(H)
% imagesc(image1)
% title('原始参考图像')
% colormap(gray)
% 
% % subplot(1,3,3);
% h1 = get(gca,'Position');
% x = h1(1);y = h1(2);w = h1(3);h  = h1(4);  %0.13,0.11,0.775,0.8015
% X1 = wzx1-rfield+X;Y1 = wzy1-rfield+Y;   %128.5,72.5
% H.Position = [672 258 260 260];
% figure(H)
% 
% axes('Position',[w/mr*(X1-(mf-1)/2)+x,h/nr*(Y1-(nf-1)/2)+y,w/mr*mf,h/nr*nf]);
% imagesc(image2)
% title('配准图像')
% colormap (gray)