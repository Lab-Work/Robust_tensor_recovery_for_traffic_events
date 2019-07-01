%this code detects the efficiency of rmc with different combinations of
%rank and observation ratio. Corruption ratio 0.1

clear;
rng('default');
rng(1);
addpath ..;
addpath PROPACK/tensor_toolbox-master ;
addpath PROPACK;

%set parameters
ni = 10; % # of run times
n1 =48; % # of trials for observation ratio;
I = 70; % size of rensor is (I,I,I)
lambda = 1/sqrt(I); % optimazation parameter
Ratio_D =  linspace(0.4,1,n1);  % observation ratio from 30% to 100%
Ranks =  2:35;  % rank (2,2,2) to (35,35,35);
n2 = length(Ranks); % # of trials for rank;
ratio_s = 0.1; % corruption ratio

% pre-allocate
err_L = zeros(ni,n1,n2); 
Precision = zeros(ni,n1,n2);  %
Recall = zeros(ni,n1,n2);  % 


tic; % start the timer

for ii = 1:ni
    disp(['iteration= ' num2str(ii)]);
for i = 1:n2 % rank
    disp(['rank= ' num2str(Ranks(i))]);
    % simulate full observsation
    c = Ranks(i);
    [D0 ,L ,S0]= simulate_commplete(I,I,I,c,c,c,ratio_s); 
    
    parfor j = i:n1 % observation ratio
        Sigma_ = rand(size(D0));
        Sigma = Sigma_< Ratio_D(j); % keep data by ratio
        Sigma = tensor(Sigma);

        S = S0.* Sigma;
        D = D0 .* Sigma;
        Sigma_bar = Sigma_ > Ratio_D(j);
        
        [Lhat ,Shat, ~,~] = inexact_alm_rmc21(D,Sigma_bar,lambda, 1e-7, 1000);
        err_L(ii,j,i) = norm(Lhat - L)/ norm(L);

        % precision
        loc = S.data~=0;
        loc_hat = Shat.data~=0;
%         loc_hat = Shat.data >= 1e-6;

        tp = sum(sum(sum((loc==1) & (loc_hat==1))));
        fn = sum(sum(sum((loc==1) & (loc_hat==0))));
        fp = sum(sum(sum((loc==0) & (loc_hat==1))));

        Precision(ii,j,i) =  tp/(tp+fp);
        Recall(ii,j,i) = tp/(tp+fn);

        if mod( j, 10) == 0
            disp(['ratio_d: ' num2str( Ratio_D(j)) '  err of L :' num2str(err_L(ii,j,i))]);
        end
    end
 
end
end

totalT = toc;
fprintf('%d minutes and %f seconds\n', floor(totalT/60), rem(totalT,60));
