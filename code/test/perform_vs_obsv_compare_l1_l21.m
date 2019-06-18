% this code draws curve of precision as funcion of observation ratio
% one can define different corruption ratio and tensor ranks.

%%  initialize

clear;
rng('default');
rng(1);

%set parameters
ni = 10; % run ni times
n =40;  % # of trials
I = 70; % size of rensor is (I,I,I)
c = int8(2); % rank is (c,c,c)
ratio_S =  linspace(0.01,0.8,n);  % range of corruption ratio
ratio_D = 0.3; % observation
lambda1 = 1/sqrt(I); % optimazation parameter for l1
lambda2 = 1/I/0.03; % optimazation parameter for l21


% pre-allocate for l1
err_L1 = zeros(ni,n); 
Precision1 = zeros(ni,n);  %
Recall1 = zeros(ni,n);  % 
Iter1 = zeros(ni,n);  % 

% pre-allocate for l2
err_L2 = zeros(ni,n); 
Precision2 = zeros(ni,n);  %
Recall2 = zeros(ni,n);  % 
Iter2 = zeros(ni,n);  % 



tic; % start the timer

for i = 1:ni % run ni times and average
    disp(['iteration # ' num2str(i)]);
    parfor j = 1:n
        % simulate observation data
        ratio_s = ratio_S(j);
        [D, L, S, O]= simulate_rmc(I,I,I,c,c,c,ratio_D,ratio_s);
        
        % l1 decompsition
        [Lhat ,Shat, ~,iter] = inexact_alm_rmc3D(D,lambda1, 1e-7, 1000);
        err= norm(Lhat - L)/ norm(L);
%         disp(['err of L :' num2str(err)]);
        
        err_L1(i,j) = err;
        Iter1(i,j) = iter;
        % precision
        loc = S.data~=0;
        loc_hat = Shat.data~=0;
%         loc_hat = Shat.data >= 1e-6;

        tp = sum(sum(sum((loc==1) & (loc_hat==1))));
        fn = sum(sum(sum((loc==1) & (loc_hat==0))));
        fp = sum(sum(sum((loc==0) & (loc_hat==1))));

        Precision1(i,j) =  tp/(tp+fp);
        Recall1(i,j) = tp/(tp+fn);
        
                
        % l21 decompsition
        [Lhat ,Shat, ~,iter] = inexact_alm_rmc21(D,lambda2, 1e-7, 1000);
        err= norm(Lhat - L)/ norm(L);
%         disp(['err of L :' num2str(err)]);
        
        err_L2(i,j) = err;
        Iter2(i,j) = iter;
        % precision
        loc = S.data~=0;
        loc_hat = Shat.data~=0;
%         loc_hat = Shat.data >= 1e-6;

        tp = sum(sum(sum((loc==1) & (loc_hat==1))));
        fn = sum(sum(sum((loc==1) & (loc_hat==0))));
        fp = sum(sum(sum((loc==0) & (loc_hat==1))));

        Precision2(i,j) =  tp/(tp+fp);
        Recall2(i,j) = tp/(tp+fn);

        if mod( j, 10) == 1
            disp(['ratio_s: ' num2str(ratio_S(j))]);
            disp(['l1: err of L :' num2str(err_L1(i,j)) '  precision:' num2str(Precision1(i,j))  ]);
            disp(['l2: err of L :' num2str(err_L2(i,j)) '  precision:' num2str(Precision2(i,j))  ]);
        end
    end
end

totalT = toc;

err_L_avg1 = mean(err_L1);
Iter_avg1 = mean(Iter1);
Recall_avg1 = mean(Recall1);
Precision_avg1 = mean(Precision1);


err_L_avg2 = mean(err_L2);
Iter_avg2 = mean(Iter2);
Recall_avg2 = mean(Recall2);
Precision_avg2 = mean(Precision2);
%% save
save ../SimulationData/compare_0.3.mat
%% plot
 
figure();
plot(ratio_S, err_L_avg1);
hold on;
plot(ratio_S, err_L_avg2);
hold off;
xlabel('ratio_s');
ylabel('err of L');
ylim([0,1.1]);
legend('l1','l21')

figure();
plot(ratio_S, Precision_avg1);
hold on;
plot(ratio_S,  Precision_avg2);
hold off;
xlabel('ratio_s');
ylabel('Precision');
legend('l1','l21')
ylim([0,1.1]);

figure();
plot(ratio_S, Recall_avg1);
hold on;
plot(ratio_S, Recall_avg2);
hold off;
xlabel('ratio_s');
ylabel('recall');
legend('l1','l21')
ylim([0,1.1]);








