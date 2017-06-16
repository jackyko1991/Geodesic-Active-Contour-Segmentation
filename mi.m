 function MuInf = mi(image1,image2)
%   MI    互信息测度（Mutual Information）计算函数
%   method参数选择计算互信息还是归一化互信息（Normalized Mutual Information）
%   作者：lskyp    2009.07.23
%   版本：V1.0
% 联合直方图归一化
% image1 = imread('c:\matlab\drlse\1.bmp');
% image2 = imread('c:\matlab\drlse\8.bmp');
% method = '互信息';
Joint_h = Jointh_My(image1,image2);
[r,c] = size(Joint_h);    %256*256
N_Jh = Joint_h./(r*c);
Marg_A = sum(N_Jh);
Marg_B = sum(N_Jh,2);
H_A = 0;H_B = 0;
% 熵函数
for k = 1:r
    if Marg_A(k) ~= 0
        H_A = H_A+(-Marg_A(k)*log2(Marg_A(k)));
    end
end
for k = 1:c
    if Marg_B(k) ~= 0
        H_B = H_B+(-Marg_B(k)*log2(Marg_B(k)));
    end
end
H_AB = sum(sum(-N_Jh.*log2(N_Jh+(N_Jh == 0))));
 
% 输出
% if strmatch(method,'互信息')
     MuInf = H_A+H_B-H_AB;
% else
%     MuInf = (H_A+H_B)/H_AB;
% end
% disp(H_A);    %4.227
% disp(H_B);    %4.4649
% disp(H_AB);   %7.3196
% disp(MuInf);  %1.3726
