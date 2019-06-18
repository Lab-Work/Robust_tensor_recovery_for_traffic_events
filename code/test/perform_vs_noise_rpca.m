% this code draws curve of precision as funcion of corruption ratio for
% full-observation cases (rpca-21).

%%  initialize
clear;
rng('default');
rng(1);
addpath PROPACK;
addpath PROPACK/tensor_toolbox-master


%set parameters
ni = 10; % run ni times
n = 50;  % # of trials
I = 70; % size of rensor is (I,I,I)
c = 5; % rank
lambda = 1/I/0.03; % optimazation parameter for l2,1 norm
% lambda1 = 1/sqrt(I); % optimization parameter for l1 norm
lambda1 = 1/sqrt(I)*0.4;
Ratio_s =  linspace(0.01,0.6,n);  % corruption ratio from 1% to 60%


% pre-allocate
err_L = zeros(ni,n); 
Precision = zeros(ni,n);  %
Recall = zeros(ni,n);  % 
Iter = zeros(ni,n);  % 


tic; % start the timer

parfor i = 1:ni % run ni times and average   
    disp(['iteration # ' num2str(i)]);
    for j = 1:n
%         corruption ratio = Ratio_s(j);
        ratio_s = Ratio_s(j);
        [D ,L ,S]= simulate_complete(I,I,I,c,c,c,ratio_s); 
%         [Lhat ,Shat, iter] = inexact_alm_rpca21(D,lambda, 1e-7, 1000);
        [Lhat ,Shat, iter] = inexact_alm_rpca3D(D,lambda1, 1e-7, 1000);
        err_L(i,j) = norm(Lhat - L)/ norm(L);
        Iter(i,j) =iter;
        
        % precision
        loc = S.data~=0;
        loc_hat = Shat.data >= 1e-6;

        tp = sum(sum(sum((loc==1) & (loc_hat==1))));
        fn = sum(sum(sum((loc==1) & (loc_hat==0))));
        fp = sum(sum(sum((loc==0) & (loc_hat==1))));

        Precision(i,j) = tp/(tp+fp);
        Recall(i,j) = tp/(tp+fn);


        if mod( j, 10) == 1
            disp(['corruption ratio: ' num2str( Ratio_s(j)) '  err of L :' num2str(err_L(i,j)) ...
                '  precision: '  num2str(Precision(i,j))]);
        end
    end
end
totalT = toc;
fprintf('%d minutes and %f seconds\n', floor(totalT/60), rem(totalT,60));


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
save ../SimulationData/sim_rpca3D.mat
%% plot

figure(1);
plot(Ratio_s, err_L_avg,'LineWidth', 2 );
xlabel('corruption ratio');
ylabel('RE');
title({'RE v.s. corruption ratio'});
set(gca, 'fontsize', 18)
print('../Figure/err_noise3D','-depsc')
% ylim([0,1]);

figure(2);
plot(Ratio_s, Precision_avg,'LineWidth', 2 );
xlabel('corruption ratio');
ylabel('Precision');
title('Precision v.s. corruption ratio');
set(gca, 'fontsize', 18)
ylim([0,1.1]);
print('../Figure/precision_noise3D','-depsc')


figure(3);
plot(Ratio_s, Recall_avg,'LineWidth', 2 );
xlabel('corruption ratio');
ylabel('Recall');
title('Recall v.s. corruption ratio');
set(gca, 'fontsize', 18)
ylim([0,1.1]);
print('../Figure/recall_noise3D','-depsc')


figure(4);
plot(Ratio_s, Iter_avg,'LineWidth', 2 );
xlabel('corruption ratio');
ylabel('# iteration');
title('# iteration v.s. corruption ratio');
set(gca, 'fontsize', 18)
ylim([13,37]);
print('../Figure/iter_noise3D','-depsc')











