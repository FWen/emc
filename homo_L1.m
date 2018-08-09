function [Inliers] = homo_L1(data,eps)
kap = 4;
m = size(data.x1,2);
n = 8;

eps = eps*data.T2(1,1);  
        
[A0, b0, c0, d0] = genMatrixHomography(data.x1, data.x2);
[A, b] = genLinearMatrixFromQuasiconvex(A0, b0, c0, d0, eps);

C = [b; sparse(kap*m,1)]; 
B = [sparse(n,1); ones(kap*m,1)]; 
A = [A, -speye(kap*m); 
    sparse(kap*m,n), -speye(kap*m)];

K.l = size(A,1);
pars.eps = 1e-8;
pars.maxiter = 1e3;
pars.fid = 0;

[~,Y,~] = sedumi(A,-B,C,K,pars);
slack = sum(reshape(Y(n+1:end),kap,m));
Inliers = find(slack<1e-7);
