function [x,s,out] = adm_alg(A,b,kap,para);

[m,n] = size(A);
ns = m/kap;

%Initialize
if(isfield(para, 'MAX_ITER'))
    MAX_ITER = para.MAX_ITER;
else
    MAX_ITER = 1e3;
end

if(isfield(para, 'TOL'))
    ABSTOL = para.TOL;
else
    ABSTOL = 1e-10;
end

if(isfield(para, 'x0'))
    x = para.x0;
else
    x = randn(n,1);
end

if(isfield(para, 's0'))
    s = para.s0;
else
    s = zeros(ns,1);
end

if(isfield(para, 'rho'))
    rho = para.rho;
else
    rho = 1;
end

e1 = ones(kap,1);
v  = min(A*x - b - kron(s,e1), 0);
u  = zeros(m,1);
AA = inv(A'*A);
out.ex = [];out.es = [];
num_stop = 0;

for iter = 1:MAX_ITER
    xm1 = x; 
    sm1 = s;
    
    % for acceleration of the algorithm
    rho = rho * 1.04;

    % x-update using the conjugate gradient method
    x = AA*(A'*(b + kron(s,e1) + v - u/rho));
    
    % s-update
    cs = A*x - b - v + u/rho;
    s = sum(reshape(cs, kap, m/kap))'/kap;
    s(s<sqrt(2/kap/rho)) = 0;
    
    % v-update
    v = min(A*x - b - kron(s,e1) + u/rho, 0);
    
    % w-update
    u = u + rho*(A*x - b - kron(s,e1) - v);

    % Iteration gap
%     out.ex = [out.ex, norm(x-xm1)/sqrt(n)];
%     out.es = [out.es, norm(s-sm1)/sqrt(ns)];
    
    % Check for convergence       
    if (norm(s-sm1)< sqrt(ns)*ABSTOL) & (norm(x-xm1)< sqrt(n)*ABSTOL)
        num_stop = num_stop + 1;
        if num_stop==3
           break;
        end
    else
        num_stop = 0;
    end

end
