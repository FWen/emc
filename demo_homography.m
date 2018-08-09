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
EPConfig.alphaMax = 1e9;

inlier_th = 4  ; % Inlier threshold in pixel
    
imgs={'christ_church','university','bodleian','magdalen','radcliffe',...
      'AerialviewsI','Corridor','kapel', 'MertonCollegeII','MertonCollegeIII','ValbonneChurch',...
      'boat','bark','bikes','graf','trees',...
      'building4','building5','building22','building24','building28',...
       'building37','building59','building67','building199'};
   
for k=1:length(imgs)
    load(['dataset/', imgs{k}, '.mat']);
    Ns(k) = size(data.x1,2);
%     figure(1); imagesc([data.matches.im1,data.matches.im2]); axis image off;

    %-RANSAC------------------------------------------
    [ransacTheta, ransacInliers, ransacRuntime] = homographyFit(data, inlier_th, 'RANSAC', rand(8,1));
    mc(k,1) = length(ransacInliers);
    rt(k,1) = ransacRuntime;


    % L1 ---------------------------------------------
    tic; 
    [L1Inliers] = homo_L1(data,inlier_th);
    mc(k,2) = length(L1Inliers);
    rt(k,2) = toc;


    % L1-(reduced dimension) -------------------------
    tic; 
    [L1rdInliers] = homo_L1rd(data,inlier_th);
    mc(k,3) = length(L1rdInliers);
    rt(k,3) = toc;


    % Iteratively-Reweighted -------------------------
    tic; 
    [L1irwInliers] = homo_L1rd_irw(data,inlier_th,20);
    mc(k,4) = length(L1irwInliers);
    rt(k,4) = toc;


    % EP using RANSAC----------------------------------
    [eprsTheta, eprsInliers, eprsRuntime] = homographyFit(data, inlier_th, 'EP', ransacTheta, EPConfig);
    mc(k,5) = length(eprsInliers);
    rt(k,5) = eprsRuntime+ransacRuntime;


    % ADMM (Algorithm 1) initialized by L1-------------------------------
    tic; 
    [admInliers] = homo_adm(data,inlier_th);
    mc(k,6) = length(admInliers);
    rt(k,6) = toc;


    % BCD (Algorithm 2)--------------------------------
    tic; 
    [bcdInliers] = homo_bcd(data,inlier_th);
    mc(k,7) = length(bcdInliers);
    rt(k,7) = toc;

    
    %ADMM (Algorithm 1) initialized by RANSAC-----------------------------
    tic; 
    [admrsInliers] = homo_adm(data,inlier_th,ransacTheta);
    mc(k,8) = length(admrsInliers);
    rt(k,8) = toc + ransacRuntime;
    
    % BCD (Algorithm 2) initialized by RANSAC-------------------------------
    tic; 
    [bcdrsInliers] = homo_bcd(data,inlier_th,ransacTheta);
    mc(k,9) = length(bcdrsInliers);
    rt(k,9) = toc + ransacRuntime;
        
    figure;
    plot_match(data.matches, [data.matches.X1; data.matches.X2], bcdrsInliers, 0, 100);
end

mc
rt

