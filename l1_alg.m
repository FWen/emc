function [slack] = l1_alg(x,y,eps)

[m, n] = size(x);
A = zeros(2*m,n);
A(1:2:end) = x;
A(2:2:end) = -x;

b = zeros(2*m,1);
b(1:2:end) = y + eps;
b(2:2:end) = -y + eps;

C = [b; sparse(size(A,1),1)]; 

B = [sparse(n,1); ones(2*m,1)]; 
A = [A, -speye(2*m); 
    sparse(2*m,n), -speye(2*m)];

K.l = size(A,1);
pars.eps = 1e-10;
pars.maxiter = 1e3;
pars.fid = 0;

[~,Y,~] = sedumi(A,-B,C,K,pars);
slack = sum(reshape(Y(n+1:end),2,m));
slack(slack<1e-7) = 0;
