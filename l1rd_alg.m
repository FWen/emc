function [slack] = l1rd_alg(x,y,eps)

[m, n] = size(x);
A = zeros(2*m,n);
A(1:2:end) = x;
A(2:2:end) = -x;

b = zeros(2*m,1);
b(1:2:end) = y + eps;
b(2:2:end) = -y + eps;

C = [b; zeros(m,1)]; 
B = [sparse(n,1); ones(m,1)]; 
J = kron(speye(m),ones(2,1));
A = [A, -J; sparse(m,n), -speye(m)];

K.l = size(A,1);
pars.eps = 1e-10;
pars.maxiter = 1e3;
pars.fid = 0;

[~,Y,~] = sedumi(A,-B,C,K,pars);
slack = Y(n+1:end);
slack(slack<1e-7) = 0;
