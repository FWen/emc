function [x,s,out] = bcd_alg(A,b,kap,para);

[m,n] = size(A);
ns = m/kap;

% Initialize
if(isfield(para, 'MAX_ITER'))
    MAX_ITER = para.MAX_ITER;
else
    MAX_ITER = 5e3;
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

MIN_ITER = 1e3;
phi = 1e-2;
e1 = ones(kap,1);
v = min(A*x - b - kron(s,e1),0);

out.ex = [];out.es = [];
AA = inv(A'*A);
num_stop = 0;

for iter = 1:MAX_ITER
    xm1 = x; 
    sm1 = s;
	
    % for acceleration of the algorithm
     if rho<5e5
        rho = rho * 1.01;
     end
    
    % x-update
    x = AA*(A'*(b + kron(s,e1) + v));
    
    % s-update
    cs = A*x - b - v;
    s = (rho*sum(reshape(cs, kap, m/kap))'+phi*sm1)/(kap*rho+phi);
    s(s<sqrt(2/(kap*rho+phi))) = 0; 
   
    % v-update
    v = min((rho*(A*x - b - kron(s,e1)) + phi*v)/(rho+phi),0);

    % iteration gap
%     out.ex = [out.ex, norm(x-xm1)/sqrt(n)];
%     out.es = [out.es, norm(s-sm1)/sqrt(ns)];
    
    %Check for convergence
    if (iter>MIN_ITER & norm(s-sm1)< sqrt(ns)*ABSTOL) & (norm(x-xm1)< sqrt(n)*ABSTOL)
        num_stop = num_stop + 1;
        if num_stop==3
            break;
        end
    else
        num_stop = 0;
    end
    
end
