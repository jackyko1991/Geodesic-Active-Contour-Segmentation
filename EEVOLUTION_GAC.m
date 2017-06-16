function u = EEVOLUTION_GAC(initLSF, g, alf, delt, numIter)
%  EVOLUTION_GAC(initLSF, g, alf, delt, numIter)
%  EVOLUTION_GAC evolves the level set function according to geodesic active contour model.
%  initLSF:  initial level set function.
%  g: edge indicator function
%  [gx, gy]: the gradient of g;
%  alf: weight of balloon/pressure force
%  delt: time step of iteration
%  numIter: number of iterations


u=initLSF;
L = length(g);
[grow gcol] = size(g{1});
% gm = ones(size(g{1}));    %取最小值
gm = ones(size(g{1}));   
S = grow * gcol;
% 取最大值/最小值/平均值
for i=1:L
   for j = 1:S
%        
% %          gm(j) = gm(j)+ 1/8*g{i}(j);    
        gm(j) = min(gm(j),g{i}(j));    %取最小值
%         gm(j) = min(gm(j),g{i}(j));    %取最小值
   end
end
% for i = 1:S
% %     gm(i) =1/4*( min(g{1}(i),g{8}(i)) + min(g{2}(i),g{7}(i)+...
% %             min(g{3}(i),g{6}(i)+min(g{4}(i),g{5}(i))))); 
% %       gm(i) =1/2*( min(g{1}(i),g{2}(i)) + min(g{3}(i),g{4}(i)+...
% %               min(g{5}(i),g{6}(i)+min(g{7}(i),g{8}(i))))); 
%        gm(i) =1/4*( max(g{1}(i),g{2}(i)) + max(g{3}(i),g{4}(i)+...
%                max(g{5}(i),g{6}(i)+max(g{7}(i),g{8}(i))))); 
% %         gm(i) =1/2*(max(max(g{1}(i),g{3}(i)),max(g{2}(i),g{4}(i)))+...
% %                max(max(g{4}(i),g{5}(i)),max(g{6}(i),g{7}(i))));  
% end
%  

    [gx,gy]=gradient(gm);

for k=1:numIter      %numIter=10 

    K=CURVATURE(u);
    dx_f=Dx_forward(u);
    dy_f=Dy_forward(u);
    dx_b=Dx_backward(u);
    dy_b=Dy_backward(u);
    [dx_c,dy_c]=gradient(u);

    norm_grad_p_u=gradientNorm_upwind(u,'+');
    norm_grad_n_u=gradientNorm_upwind(u,'-');    
   
    u=u+delt*(gm.*K.*sqrt(dx_c.^2+dy_c.^2) ...
        + alf*(max(gm,0).*norm_grad_p_u + min(gm,0).*norm_grad_n_u) ...
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
