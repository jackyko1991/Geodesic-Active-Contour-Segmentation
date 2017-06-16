function fx=Dx_forward(f)
% Dx_forward(f) computes forward difference with respect to x (or j, the
% idex of column) for coordinate system as below
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
fx=zeros(nr,nc);
fx(:,1:nc-1)=f(:,2:nc)-f(:,1:nc-1);
fx(:,nc) = f(:,nc)-f(:,nc-1);  % for the last column, use the backward difference as the forward difference


