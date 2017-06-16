function d = gradientNorm_upwind(u,sign)
%   gradientNorm_upwind(u,sign)
%   Chunming Li  04/26/2004
%   Copyright (c) 2001-2005 by Chunming Li

dx_f=Dx_forward(u);
dy_f=Dy_forward(u);
dx_b=Dx_backward(u);
dy_b=Dy_backward(u);

if sign=='+'
    d=sqrt((max(dx_b,0)).^2+(min(dx_f,0)).^2 ...
        + (max(dy_b,0)).^2+(min(dy_f,0)).^2);
else
    d=sqrt((max(dx_f,0)).^2+(min(dx_b,0)).^2 ...
        + (max(dy_f,0)).^2+(min(dy_b,0)).^2); 
end