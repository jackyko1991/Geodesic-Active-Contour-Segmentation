function u = EVOLUTION_GAC(initLSF, g, alf, delt, numIter)
%  EVOLUTION_GAC(initLSF, g, alf, delt, numIter)
%  EVOLUTION_GAC evolves the level set function according to geodesic active contour model.
%  initLSF:  initial level set function.
%  g: edge indicator function
%  [gx, gy]: the gradient of g;
%  alf: weight of balloon/pressure force
%  delt: time step of iteration
%  numIter: number of iterations
%
% Author: Chunming Li
% email: li_chunming@hotmail.com

u=initLSF;
[gx,gy]=gradient(g);
for k=1:numIter      %numIter=10 

    K=CURVATURE(u);
    dx_f=Dx_forward(u);
    dy_f=Dy_forward(u);
    dx_b=Dx_backward(u);
    dy_b=Dy_backward(u);
    [dx_c,dy_c]=gradient(u);

    norm_grad_p_u=gradientNorm_upwind(u,'+');
    norm_grad_n_u=gradientNorm_upwind(u,'-');    
   
    u=u+delt*(g.*K.*sqrt(dx_c.^2+dy_c.^2) ...
        + alf*(max(g,0).*norm_grad_p_u + min(g,0).*norm_grad_n_u) ...
        + (max(gx,0).*dx_b + min(gx,0).*dx_f ...
        + max(gy,0).*dy_b + min(gy,0).*dy_f));    %image energe and evolution style 
                                                %该公式里面的g会小于0吗？
end


function K = CURVATURE(f)   %？？？？？why forward -> backward???
epsilon=1e-10;
fx_f = Dx_forward(f);
fy_f = Dy_forward(f);
ax = fx_f./sqrt(fx_f.^2+ fy_f.^2+epsilon);
ay = fy_f./sqrt(fx_f.^2 + fy_f.^2 + epsilon);
axx = Dx_backward(ax);
ayy = Dy_backward(ay);
K = axx + ayy;
