function [slack] = max_con_adm(x,y,epsilon,x0)

kap = 2;
[m, n] = size(x);
A1 = zeros(2*m,n);
A1(1:2:end) = x;
A1(2:2:end) = -x;

b = zeros(2*m,1);
b(1:2:end) = y + epsilon;
b(2:2:end) = -y + epsilon;


%Initialization
if nargin>=4
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
para.MAX_ITER = 2e3;
para.TOL = 1e-12;
para.rho=1;
[~,slack,~] = adm_alg(A1,b,2,para);
