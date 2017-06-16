function Trad = IsInTriangleArea(p,p1,p2,p3)
% 判断点P是否在三角形区域内
% 用三角形面积来判断。原三个顶点p1,p2,p3,新点p.
% (s(p,p1,p2)+s(p,p1,p3)+s(p,p3,p2))=s(p1,p2,p3) 
s0 = s(p1,p2,p3);
s1 = s(p,p1,p2);
s2 = s(p,p1,p3);
s3 = s(p,p2,p3);
if round(s0*10) == round((s1+s2+s3)*10)
    Trad = 1;
else 
    Trad = 0;
end