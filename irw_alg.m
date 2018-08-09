function [slack] = irw_alg(x,y,eps,iter)

[m, n] = size(x);
A1 = zeros(2*m,n);
A1(1:2:end) = x;
A1(2:2:end) = -x;

b = zeros(2*m,1);
b(1:2:end) = y + eps;
b(2:2:end) = -y + eps;

C = [b; sparse(m,1)]; 
B = [sparse(n,1); ones(m,1)]; 
J = kron(speye(m),ones(2,1));
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
slack(slack<1e-10)=0;
