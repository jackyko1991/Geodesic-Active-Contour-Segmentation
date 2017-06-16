function phi = VEVOLUTION_GAC(initLSF, gm, alfa, delt, numIter)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
%  EVOLUTION_GAC evolves the level set function according to geodesic active contour model.
%  initLSF:  initial level set function.
%  g: edge indicator function
%  [gx, gy]: the gradient of g;
%  alf: weight of balloon/pressure force
%  delt: time step of iteration
%  numIter: number of iterations


% u=initLSF;
phi = initLSF;
[gx,gy]=gradient(gm);


for k=1:numIter      %numIter=10 
    [phi_x,phi_y]=gradient(phi);
    s=sqrt(phi_x.^2 + phi_y.^2);
    smallNumber=1e-10;  
    Nx=phi_x./(s+smallNumber); % add a small positive number to avoid division by zero
    Ny=phi_y./(s+smallNumber);
    curvature=div(Nx,Ny);
    diracPhi=Dirac(phi,1.5);   %结果与除以s有着很大的关系
    areaTerm=diracPhi.*gm; % balloon/pressure force
    edgeTerm=diracPhi.*(gx.*Nx+gy.*Ny) + diracPhi.*gm.*curvature;  
%     distRegTerm = 4*del2(phi)-curvature; %常规项
%       distRegTerm=distReg_p2(phi);  
    phi=phi + delt*(5*edgeTerm + alfa*areaTerm);
    
end


function f = div(nx,ny)
[nxx,junk]=gradient(nx);  
[junk,nyy]=gradient(ny);
f=nxx+nyy;

function f = Dirac(x, sigma)
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;

function f = distReg_p2(phi)
% compute the distance regularization term with the double-well potential p2 in eqaution (16)
[phi_x,phi_y]=gradient(phi);
s=sqrt(phi_x.^2 + phi_y.^2);
a=(s>=0) & (s<=1);
b=(s>1);
ps=a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  % compute first order derivative of the double-well potential p2 in eqaution (16)
dps=((ps~=0).*ps+(ps==0))./((s~=0).*s+(s==0));  % compute d_p(s)=p'(s)/s in equation (10). As s-->0, we have d_p(s)-->1 according to equation (18)
f = div(dps.*phi_x - phi_x, dps.*phi_y - phi_y) + 4*del2(phi);  
