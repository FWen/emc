function [Inliers] = homo_L1rd_irw(data,eps,iter)

kap = 4;
m = size(data.x1,2);
n = 8;

% x1 = data.x1;
% x2 = data.x2;        
% T2 = data.T2;
eps = eps*data.T2(1,1);          
[A0, b0, c0, d0] = genMatrixHomography(data.x1, data.x2);
[A1, b] = genLinearMatrixFromQuasiconvex(A0, b0, c0, d0, eps);    

C = [b; zeros(m,1)]; 
B = [sparse(n,1); ones(m,1)]; 
J = kron(speye(m),ones(kap,1));
A = [A1, -J; sparse(m,n), -speye(m)];

K.l = size(A,1);
pars.eps = 1e-8;
pars.maxiter = 1e3;
pars.fid = 0;

[~,Y,~] = sedumi(A,-B,C,K,pars);
slack = Y(n+1:end);

for k=1:iter-1
    w = (max(slack,0)+1e-3).^(-0.9);
    B = [sparse(n,1); w];
    [~,Y,~] = sedumi(A,-B,C,K,pars);
    slack = Y(n+1:end);
end
Inliers = find(slack<1e-10);
