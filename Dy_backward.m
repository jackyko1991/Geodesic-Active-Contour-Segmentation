function fy=Dy_backward(f);
% Dy_backward(f) computes backward difference with respect to y (or i, the
% idex of row) for coordinate system as below
% 
%
%           x(j)
%        |--------->
%    y(i)|
%        |
%        |
%        V
%
% Copyright (c) 2001--2006 by Chunming Li
% Author: Chunming Li, 11/02/2003
% Revision by Chunming Li 7/24/2005

[nr,nc]=size(f);
fy=zeros(nr,nc);
fy(2:nr,:)=f(2:nr,:)-f(1:nr-1,:);

fy(1,:)=f(2,:)-f(1,:); % for the first row, use the forward difference as the backward difference
