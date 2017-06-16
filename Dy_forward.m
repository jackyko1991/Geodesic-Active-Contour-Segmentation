function fy=Dy_forward(f)
% Dy_forward(f) computes forward difference with respect to y (or i, the
% idex of row) for coordinate system as below
% 
%
%           x(j)
%        |---------->
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
 fy(1:nr-1,:)=f(2:nr,:)-f(1:nr-1,:);
fy(nr,:)=f(nr,:)-f(nr-1,:);  % for the last row, use the backward difference as the forward difference

