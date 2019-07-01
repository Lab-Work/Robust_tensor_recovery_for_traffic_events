% this code draws curve of precision as funcion of observation ratio
% one can define different corruption ratio and tensor ranks.
% use this for rank-plot and noise-plot

%%  initialize

addpath ..;
addpath PROPACK/tensor_toolbox-master
addpath PROPACK;

clear;
rng('default');
rng(1);

%set parameters
ni = 1; % run ni times
n =10;  % # of trials
I = 70; % size of rensor is (I,I,I)
c = int8(5); % rank is (c,c,c)
Ratio_D =  linspace(0.1,1,n);  % range of observation ratio
ratio_s = 0.1; % corruption ratio
lambda = 1/I/0.03; % optimazation parameter
lambda1 = 1/sqrt(I)*0.4;


% pre-allocate
err_L = zeros(ni,n); 
Precision = zeros(ni,n);  %
Recall = zeros(ni,n);  % 
Iter = zeros(ni,n);  % 


tic; % start the timer

for i = 1:ni % run ni times and average
    disp(['iteration # ' num2str(i)]);
    % simulate full observsation
    [D0 ,L ,S0]= simulate_complete(I,I,I,c,c,c,ratio_s); 
    parfor j = 1:n
        % partial observation
        Sigma_ = rand(size(D0));
        Sigma = Sigma_< Ratio_D(j); % keep data by ratio
        Sigma_bar = Sigma_ >= Ratio_D(j); % index for missing data
        Sigma = tensor(Sigma);

        S = S0.* Sigma;
        D = D0 .* Sigma;
%         [Lhat ,Shat, ~,iter] = inexact_alm_rmc21_v3(D,Sigma_bar,lambda, 1e-7, 1000);
        [Lhat ,Shat, ~,iter] = inexact_alm_rmc3D(D,Sigma_bar,lambda1, 1e-7, 1000);
        err= norm(Lhat - L)/ norm(L);
%         disp(['err of L :' num2str(err)]);
        
        err_L(i,j) = err;
        Iter(i,j) = iter;
        % precision
        loc = S.data~=0;
        loc_hat = Shat.data~=0;
%         loc_hat = Shat.data >= 1e-6;

        tp = sum(sum(sum((loc==1) & (loc_hat==1))));
        fn = sum(sum(sum((loc==1) & (loc_hat==0))));
        fp = sum(sum(sum((loc==0) & (loc_hat==1))));

        Precision(i,j) =  tp/(tp+fp);
        Recall(i,j) = tp/(tp+fn);

        if mod( j, 10) == 0
            disp(['ratio_d: ' num2str( Ratio_D(j)) '  err of L :' num2str(err)]);
        end
    end
end

totalT = toc;

if ni > 1
    err_L_avg = mean(err_L);
    Iter_avg = mean(Iter);
    Recall_avg = mean(Recall);
    Precision_avg = mean(Precision); 
else
    err_L_avg = (err_L);
    Iter_avg = (Iter);
    Recall_avg = (Recall);
    Precision_avg = (Precision);
end
%% save
save ../../SimulationData/sim_rmc3D_0.05_70_5.mat
%% plot
 
% figure();
% plot(Ratio_D, err_L_avg);
% xlabel('observation ratio');
% ylabel('err of L');
% % ylim([0,1.1]);
figure()
plot(Ratio_D, err_L_avg,'LineWidth', 2 );
xlabel('observation ratio');
ylabel('RE');
title({'RE v.s. observation ratio','for different corruption ratios'});
% ylim([0,1]);
set(gca, 'fontsize', 18)

figure();
plot(Ratio_D, Precision_avg);
yyaxis right
ylabel('Precision');
hold on
plot(Ratio_D, Recall_avg);
xlabel('ratio_ s');
yyaxis left
ylabel('Recall');
hold off


figure();
plot(Ratio_D, Iter_avg);
xlabel('ratio_s');
ylabel('iteration#');







