function [Inliers] = homo_bcd(data,eps,x0)

kap = 4;
m = size(data.x1,2);
n = 8;
eps = eps*data.T2(1,1);%normlization
[A, b, c, d] = genMatrixHomography(data.x1, data.x2);
[A1, b] = genLinearMatrixFromQuasiconvex(A, b, c, d, eps);    

%Initialization
if nargin>=3
    if (length(x0)==9)
        x0 = x0./x0(end);             
        x0 = x0(1:n);
    end                   
    para.x0 = x0;
    
    s = max(reshape(A1*x0-b,kap,m)).';
    para.s0 = max(s,0);
else
   %--Obtain an initialization by the L1-RD
    C = [b; zeros(m,1)]; 
    B = [sparse(n,1); ones(m,1)]; 
    J = kron(speye(m),ones(kap,1));
    A = [A1, -J; sparse(m,n), -speye(m)];

    K.l = size(A,1);
    pars.eps = 1e-7;
    pars.maxiter = 1e3;
    pars.fid = 0;
    [~,Y,~] = sedumi(A,-B,C,K,pars);
    para.x0 = Y(1:n);
    para.s0 = Y(n+1:end);
end

%-- Proximal BCD algorithm--
para.MAX_ITER = 5e3;
para.TOL = 1e-8;
para.rho=10;
[~,slack,~] = bcd_alg(A1,b,kap,para);
Inliers = find(slack==0);
