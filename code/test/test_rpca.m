
% This code Uses simulate_complete and inexact_al_rpca21.m
% to conduct robust tensor RPCA
% can choose to plot the color map of the tensor (expanded in mode 1)
%
% Yue Hu, Jun 2019. Questions? yue.hu@vadnerbilt.edu;

%% simulate
clear;

addpath ../*;
addpath ../PROPACK;
addpath ../PROPACK/tensor_toolbox-master ;

Plot = 0;  % if Plot == 1, then plot the result
rng('default');
rng(8);

ratio_s = 0.05;

% simulate_4D for 4D and simulate_complete for 3D
% [D ,L ,S]= simulate_complete(219,168,16,100,20,2,ratio_s);
% [D ,L ,S]= simulate_4D(219,24,7,16,8,2,3,2,ratio_s);
% [D ,L ,S]= simulate_4D(27,27,27,27,8,2,3,2,ratio_s);

I = 70;
c = 0.1 * I;
[D ,L ,S]= simulate_complete(I,I,I,c,c,c,ratio_s);
%% RPCA
lambda = 1/I/0.03;
% lambda = 1/sqrt(log(I)); 
lambda1 = 1/sqrt(I)*0.4;

fprintf('lambda = %1.4f \n',lambda)

tic;
[Lhat ,Shat,iter] = inexact_alm_rpca21(D, lambda, 1e-7, 1000);

% [Lhat ,Shat,iter] = inexact_alm_rpca3D(D, lambda1, 1e-7, 1000);
time = toc;

rss3 = norm(Lhat - L)/ norm(L);
rss4 = norm(Shat - S) / norm(S);
Size = numel(double(D));
tol_spar = 1e-7;  % tolerance for non-zero in Shat;
Spar = sum(sum(sum(sum(abs(double(Shat)) >= tol_spar))))/Size;  % percentage of non-zero elements in S
disp(['Estimated Sparcity: ', num2str(Spar)])
fprintf('For L2,1 norm, residual for low rank matrix is %e\n',rss3)
fprintf('residual for sparse matrix is %e\n',rss4)
fprintf('total iteration: %d\n',iter)
fprintf('total time elapsed: %.3f\n',time)

% precision
loc = S.data~=0;
loc_hat = abs(Shat.data)  >= tol_spar ;

tp = sum(sum(sum(sum((loc==1) & (loc_hat==1)))));
fn = sum(sum(sum(sum((loc==1) & (loc_hat==0)))));
fp = sum(sum(sum(sum((loc==0) & (loc_hat==1)))));


Precision = tp/(tp+fp);
Recall = tp/(tp+fn);
disp(['Precision: ', num2str(Precision)])
disp(['Recall: ', num2str(Recall)])

fprintf('----------------\n')



%% plot
if Plot == 1

    D_f = tenmat(D,1);
    Shat_f = tenmat(Shat,1);
    Lhat_f = tenmat(Lhat,1);
    L_f = tenmat(L,1);
    S_f = tenmat(S,1);

    figure('position',[0,0,10000,200]);
    im1 = image(D_f.data,'CDataMapping','scaled');
    caxis([-0.015 0.015]);
    colorbar
    title('Observation data');


    figure('position',[0,200,10000,200]);
    im1 = image(L_f.data,'CDataMapping','scaled');
    caxis([-0.015 0.015]);
    colorbar
    title('Original Low Rank data');

    figure('position',[0,400,10000,200]);
    im3 = image(Lhat_f.data,'CDataMapping','scaled');
    caxis([-0.015 0.015]);
    colorbar;
    title('Predicted Low Rank data');

    figure('position',[0,600,10000,200]);
    im2 = image(Shat_f.data,'CDataMapping','scaled');
    caxis([-0.015 0.015]);
    colorbar;
    title('Predicted Sparsity');

    figure('position',[0,800,10000,200]);
    im3 = image(S_f.data,'CDataMapping','scaled');
    caxis([-0.015 0.015]);
    colorbar;
    title('Original Sparcity');
end




