function fx=Dx_backward(f)
% Dx_backward(f) computes backward difference with respect to x (or j, the
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
fx(:,2:nc)=f(:,2:nc)-f(:,1:nc-1);
fx(:,1)=f(:,2)-f(:,1); % for the first column, use the forward difference as the backwar

% d difference

