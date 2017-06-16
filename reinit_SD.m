function [phi] = reinit_SD(phi, dx, dy, alpha, iterations)
% 重新初始化水平集（满足符号距离函数）
% function [phi] = reinit_SD(phi, dx, dy, alpha, iterations)
%
% Reinitializes phi into a signed distance function while preserving
% the zero level set (the interface or the curve).
%
% dx and dy are the resolution of the grid at x and y dimensions.
% alpha is a constant for calculating the euler step (dt). Should
% be between 0 and 1. 0.5 is quite safe whereas 0.9 can be risky.
% iterations specifies the number of iterations before the function returns.
% These correspond to 1st, 2nd, 3rd and 5th order accurate schemes 
% for calculating the derivative of phi.
%
% Author: Baris Sumengen  sumengen@ece.ucsb.edu
% http://vision.ece.ucsb.edu/~sumengen/
%

init_normal = @init_normal_ENO1;
evolve_normal = @evolve_normal_ENO1;
	
S_phi_0 = phi./sqrt(phi.^2 + dx.^2);
% S_phi_0 = phi./sqrt(dy.^2 + dx.^2);

Vn_ext = feval(init_normal, S_phi_0);   
it=0;
t=0;
while(it < iterations)
	[delta_normal, H1_abs, H2_abs] = feval(evolve_normal, phi, dx, dy, Vn_ext);
	dt = get_dt_normal(alpha, dx, dy, H1_abs, H2_abs);
	phi = phi + dt*(S_phi_0 - delta_normal);
	it = it+1;
	t = t+dt;
end


function [dt] = get_dt_normal(alpha, dx, dy, H1_abs, H2_abs)
if alpha <= 0 | alpha >= 1 
    error('alpha needs to be between 0 and 1!');
end

maxs = max(max(H1_abs/dx + H2_abs/dy));
dt = alpha/(maxs+(maxs==0));

function [der] = select_der_normal(Vn, der_minus, der_plus)

if size(der_minus) ~= size(der_plus) | size(der_plus) ~= size(Vn)
    error('plus, minus derivative vectors and normal force (Vn) need to be of equal length!');
end

der = zeros(size(der_plus));

for i=1:numel(Vn)
	Vn_der_m = Vn(i)*der_minus(i);
	Vn_der_p = Vn(i)*der_plus(i);
	if Vn_der_m <= 0 & Vn_der_p <= 0
		der(i) = der_plus(i);
	elseif Vn_der_m >= 0 & Vn_der_p >= 0
		der(i) = der_minus(i);
	elseif Vn_der_m <= 0 & Vn_der_p >= 0
		der(i) = 0;
	elseif Vn_der_m >= 0 & Vn_der_p <= 0
		if abs(Vn_der_p) >= abs(Vn_der_m)
			der(i) = der_plus(i);
		else
			der(i) = der_minus(i);
		end
	end
end


function [data_x] = der_ENO1_minus(data, dx)
data_x = zeros(size(data));

% extrapolate the beginning and end points of data
data(1) = 2*data(2)-data(3);
data(end) = 2*data(end-1)-data(end-2);

data_x(2:end-1) = (data(2:end-1)-data(1:end-2))/dx;


function [data_x] = der_ENO1_plus(data, dx)
data_x = zeros(size(data));

% extrapolate the beginning and end points of data
data(1) = 2*data(2)-data(3);
data(end) = 2*data(end-1)-data(end-2);

data_x(2:end-1) = (data(3:end)-data(2:end-1))/dx;


function [delta, H1_abs, H2_abs] = evolve_normal_ENO1(phi, dx, dy, Vn)

delta = zeros(size(phi)+2);
data_ext = zeros(size(phi)+2);
data_ext(2:end-1,2:end-1) = phi;

% Calculate the derivatives (both + and -)
phi_x_minus = zeros(size(phi)+2);
phi_x_plus = zeros(size(phi)+2);
phi_y_minus = zeros(size(phi)+2);
phi_y_plus = zeros(size(phi)+2);
phi_x = zeros(size(phi)+2);
phi_y = zeros(size(phi)+2);
% first scan the rows
for i=1:size(phi,1)
	phi_x_minus(i+1,:) = der_ENO1_minus(data_ext(i+1,:), dx);	
	phi_x_plus(i+1,:) = der_ENO1_plus(data_ext(i+1,:), dx);	
	phi_x(i+1,:) = select_der_normal(Vn(i+1,:), phi_x_minus(i+1,:), phi_x_plus(i+1,:));
end

% then scan the columns
for j=1:size(phi,2)
	phi_y_minus(:,j+1) = der_ENO1_minus(data_ext(:,j+1), dy);	
	phi_y_plus(:,j+1) = der_ENO1_plus(data_ext(:,j+1), dy);	
	phi_y(:,j+1) = select_der_normal(Vn(:,j+1), phi_y_minus(:,j+1), phi_y_plus(:,j+1));
end

abs_grad_phi = sqrt(phi_x.^2 + phi_y.^2);

H1_abs = abs(Vn.*phi_x.^2 ./ (abs_grad_phi+dx*dx*(abs_grad_phi == 0)));
H2_abs = abs(Vn.*phi_y.^2 ./ (abs_grad_phi+dx*dx*(abs_grad_phi == 0)));
H1_abs = H1_abs(2:end-1,2:end-1);
H2_abs = H2_abs(2:end-1,2:end-1);

delta = Vn.*abs_grad_phi;
delta = delta(2:end-1,2:end-1);

function [Vn_ext] = init_normal_ENO1(Vn)
Vn_ext = zeros(size(Vn)+2);
Vn_ext(2:end-1,2:end-1) = Vn;







