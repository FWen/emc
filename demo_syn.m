clear all; clc;

addpath('src/');
% addpath('C:\gurobi801\win64\matlab\');

solver = 'sedumi'; % solver = 'gurobi' will make gurobi as LP solver    
%the EP algorithm using 'gurobi' is several times faster than uisng 'sedumi',
%but for the  EP algorithm, uisng 'sedumi' can yield better performance than using 'gurobi'
EPConfig.lpsolver = prepareSolver(solver);
EPConfig.solver = solver;

N = 500;   % Number of points
d = 8;     % data dimension
sig = 0.1; % inlier variance used to generate Gaussian noise
osig = 1;  % outlier vanrance
balanceData = 1; % If balanceData = 1, outliers are evenly distributed on both sides of the hyperplane.
                 % Otherwise, if balanaceData = 0, outliers are distributed only on positive side of the hyperplane.

inlier_th = 0.1;   % Inlier threshold

EPConfig.maxAlpha = 1e10; %safeguard to prevent alpha to go astray
EPConfig.QThresh = 1e-4; % As Q(z) is smaller than QThresh, quit and return the result
EPConfig.alpha = 0.5;   % Initial alpha
EPConfig.kappa = 5;     % Increment rate

out_ratio = [0,5,10:10:60];

for nr=1:length(out_ratio)
    outlierP = out_ratio(nr); % Outlier percentage
    
    for k=1:20
%         [outlierP,k];

        [A,y,xt] = genRandomLinearData(N, d, sig, osig, outlierP, balanceData);
        
        tic; % L1
        [s_l1] = l1_alg(A,y,inlier_th);
        mc(1,k) = N - sum(s_l1>0);
        runtime(1,k) = toc;

        tic;  % L1-RD
        [s_l1rd] = l1rd_alg(A,y,inlier_th);
        mc(2,k) = N - sum(s_l1rd>0);
        runtime(2,k) = toc;
        
        tic; % iterativelt reweighted method
        [s_irw] = irw_alg(A,y,inlier_th,20);
        mc(3,k) = N - sum(s_irw>0);
        runtime(3,k) = toc;
        
        tic; % ADMM
        [s_adm] = max_con_adm(A,y,inlier_th);
        mc(4,k) = N - sum(s_adm>0);
        runtime(4,k) = toc;

        tic; % BCD
        [s_bcd] = max_con_bcd(A,y,inlier_th);
        mc(5,k) = N - sum(s_bcd>0);
        runtime(5,k) = toc;
               
        
      %---RANSAC--------------------------------------------------------------------
        [ransacTheta, ransacInliers, ransacRuntime ] = linearFit(A, y, inlier_th, 'RANSAC', randn(1,d));
        mc(6,k) = ransacInliers;
        runtime(6,k) = ransacRuntime;

        %---Exact Penalty (EP) with LEAST SQUARE starting point-----------------------
        [lsqTheta, lsqInliers, lsqRuntime ] = linearFit(A, y, inlier_th, 'LSQ', randn(1,d), EPConfig);
        [~, eplsqInliers, eplsqRuntime] = linearFit(A, y, inlier_th, 'EP', lsqTheta, EPConfig);
        mc(7,k) = eplsqInliers;
        runtime(7,k) = lsqRuntime+eplsqRuntime;

        %---Exact Penalty (EP)  with RANSAC starting point-----------------------
        [~, eprsInliers, eprsRuntime] = linearFit(A, y, inlier_th, 'EP', ransacTheta, EPConfig);
        mc(8,k) = eprsInliers;
        runtime(8,k) = ransacRuntime+eprsRuntime;
        
    end
    av_mc(:,nr) = mean(mc,2)
    av_rt(:,nr) = mean(runtime,2)
end

figure(11);
plot(out_ratio,av_mc(6,:),'g--',out_ratio,av_mc(1,:),'b--o',out_ratio,av_mc(2,:),'b--+',out_ratio,av_mc(3,:),'b--*',...
     out_ratio,av_mc(7,:),'g--x',out_ratio,av_mc(8,:),'g--o',...
     out_ratio,av_mc(4,:),'r--*',out_ratio,av_mc(5,:),'r--^');
xlabel('Outlier ratio (%)'); grid on;
ylabel('Consensus Size')
legend('RANSAC','L1','L1-RD','IRWLq','EP-LSQ','EP-RS','ADMM','BCD');

figure(2);
semilogy(out_ratio,av_rt(6,:),'g--',out_ratio,av_rt(1,:),'b--o',out_ratio,av_rt(2,:),'b--+',out_ratio,av_rt(3,:),'b--*',...
     out_ratio,av_rt(7,:),'g--x',out_ratio,av_rt(8,:),'g--o',...
     out_ratio,av_rt(4,:),'r--*',out_ratio,av_rt(5,:),'r--^');
xlabel('Outlier ratio (%)'); grid on;
ylabel('Runtime (seconds)')
legend('RANSAC','L1','L1-RD','IRWLq','EP-LSQ','EP-RS','ADMM','BCD');
