clear all; close all;

addpath('src/');

solver = 'sedumi'; % solver = 'gurobi' will make gurobi as LP solver    
%the EP algorithm using 'gurobi' is several times faster than uisng 'sedumi',
%but for the  EP algorithm, uisng 'sedumi' can yield better performance than using 'gurobi'
EPConfig.lpsolver = prepareSolver(solver);
EPConfig.solver = solver;
EPConfig.QThresh = 1e-5;
EPConfig.alpha = 10; 
EPConfig.kappa = 1.5;
EPConfig.maxAlpha = 1e9;

th_pixel = 2; % Inlier threshold (pixels)
d = 6;

imgs={'boat','boat2','bark','bark2','bikes','bikes2','graf','graf2','trees','trees2','wall','wall2'};
   
for k=1:length(imgs)
    load(['aff_dataset/', imgs{k}, '.mat']);
    N = size(data.x1,2);
    Ns(k) = N;
    inlier_th = th_pixel*data.T2(1,1);
%     figure(1); imagesc([data.matches.im1,data.matches.im2]); axis image off;

    [A, y] = genMatrix_Affinity(data.x1, data.x2);
    
    %-RANSAC------------------------------------------
    [ransacTheta, ransacInliers, ransacRuntime ] = linearFit_homo(A, y, inlier_th, 'RANSAC', randn(1,d));
    mc(k,1) = length(ransacInliers);
    rt(k,1) = ransacRuntime;

    
    % L1 ---------------------------------------------   
    tic; % L1
    [s_l1] = l1_alg(A,y,inlier_th);
    L1Inliers = find(sum(reshape(s_l1,2,N))==0);
    mc(k,2) = length(L1Inliers);
    rt(k,2) = toc;


    % L1-(reduced dimension) L1-RD -------------------------
    tic;  
    [s_l1rd] = l1rd_alg(A,y,inlier_th);
    L1rdInliers = find(sum(reshape(s_l1rd,2,N))==0);
    mc(k,3) = length(L1rdInliers);
    rt(k,3) = toc;

    
    tic; % Iteratively reweighted method
    [s_irw] = irw_alg(A,y,inlier_th,20);
    L1irwInliers = find(sum(reshape(s_irw,2,N))==0);
    mc(k,4) = length(L1irwInliers);
    rt(k,4) = toc;
   

    % EP using RANSAC----------------------------------
    [~, eprsInliers, eprsRuntime] = linearFit_homo(A, y, inlier_th, 'EP', ransacTheta, EPConfig);
    mc(k,5) = length(eprsInliers);
    rt(k,5) = eprsRuntime+ransacRuntime;           

    
    % ADMM (Algorithm 1) initialized by L1-------------------------------
    tic; % ADMM
    [s_adm] = max_con_adm(A,y,inlier_th);
    admInliers = find(sum(reshape(s_adm,2,N))==0);
    mc(k,6) = length(admInliers);
    rt(k,6) = toc;

    
    tic; % BCD initialized by L1
    [s_bcd] = max_con_bcd(A,y,inlier_th);
    bcdInliers = find(sum(reshape(s_bcd,2,N))==0);
    mc(k,7) = length(bcdInliers);
    rt(k,7) = toc;
    
   
    % ADMM (Algorithm 1) initialized by RANSAC-----------------------------
    tic; 
    [s_admrs] = max_con_adm(A,y,inlier_th,ransacTheta);
    admrsInliers = find(sum(reshape(s_admrs,2,N))==0);
    mc(k,8) = length(admrsInliers);
    rt(k,8) = toc + ransacRuntime;
    
   
    % BCD (Algorithm 2) initialized by RANSAC-------------------------------
    tic; 
    [s_bcdrs] = max_con_bcd(A,y,inlier_th,ransacTheta);
    bcdrsInliers = find(sum(reshape(s_bcdrs,2,N))==0);
    mc(k,9) = length(bcdrsInliers)
    rt(k,9) = toc + ransacRuntime;
    

    figure;
    plot_match(data.matches, [data.matches.X1; data.matches.X2], bcdrsInliers, 0, 100);
end

mc
rt
