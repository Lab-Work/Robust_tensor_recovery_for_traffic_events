% This code tests tensor robust completion under coloumn-wise corruotion.


%% simulate
clear;

addpath ../*;
addpath ../PROPACK;
addpath ../PROPACK/tensor_toolbox-master ;

%set seed
rng('default');
rng(8);
I = 70; % size of rensor is (I,I,I)

c = 5; % rank is (c,c,c)
lambda = 1/I/0.03; % optimazation parameter
% lambda = 1/sqrt(log(I)); 

lambda1 = 1/sqrt(I)*0.4;

[D, L, S,Sigma_bar, O]= simulate_rmc(I,I,I,c,c,c,1,0.2);
%% rpca
% [Lhat ,Shat] = inexact_alm_rpca21(D, lambda, 1e-7, 1000);
% rss = norm(Lhat - L)/ norm(L);
% rss2 = norm(Shat - S) / norm(S);
% fprintf('For rpca_21, residual for low rank matrix is %e\n',rss)
% fprintf('residual for sparse matrix is %e\n',rss2)


%% rmc

tic;
% [Lhat2,Shat, O_hat,iter] = inexact_alm_rmc21_v3(D,Sigma_bar,lambda, 1e-7, 1000);
% 
[Lhat2,Shat, O_hat,iter] = inexact_alm_rmc3D(D, Sigma_bar, lambda1, 1e-7, 1000);
toc;
rss3 = norm(Lhat2 - L)/ norm(L);
rss4 = norm(Shat - S) / norm(S);
Size = numel(double(D));
tol_spar = 1e-7; % tolerance for Shat.
Spar = sum(sum(sum(abs(double(Shat)) >= tol_spar)))/Size;  % percentage of non-zero elements in S
disp(['Estimated Sparcity: ', num2str(Spar)])
fprintf('For rmc_21, residual for low rank matrix is %e\n',rss3)
fprintf('residual for sparse matrix is %e\n',rss4)
fprintf('total iteration is %d\n',iter)

% precision
loc = S.data~=0;
loc1 = Shat.data~=0;
loc2 = abs(Shat.data) >= tol_spar;

tp = sum(sum(sum((loc==1) & (loc2==1))));
fn = sum(sum(sum((loc==1) & (loc2==0))));
fp = sum(sum(sum((loc==0) & (loc2==1))));

precision = tp/(tp+fp);
recall = tp/(tp+fn);
disp(['Precision: ', num2str(precision)])
disp(['Recall: ', num2str(recall)])

fprintf('----------------\n')

