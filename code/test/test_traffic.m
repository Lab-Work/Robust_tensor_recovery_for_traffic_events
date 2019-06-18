% this file tests robust tensor completion for realworld traffic data
% put your data in '../traffic_data folder
addpath PROPACK/tensor_toolbox-master
addpath ..

clear;
load('../traffic_data/Obs_Davidson_W17.mat')

% deal with daylight time, insert row# 1659
Obs3 = [Obs2(1:1659,:); Obs2(1659:2855,:)];
Obs_rat = 1- sum(sum(isnan(Obs3)))/numel(Obs3);  % observation ratio
disp(['Observation ratio: ', num2str(Obs_rat)])

% construct observation matrix into tensor fromat
Obs = Obs3.' ;   % for index*column,  convert from hour*link into link*hour

nl = size(Obs,1);        % #links
nh = 7*24;               % #hours
nw = size(Obs,2)/nh;     % #weeks

D_plot = tensor(Obs,[nl nh nw]);   % nan stays the same for plot; 3D
Obs_Nan = Obs;  % nan stays the same for plot; 2D

% index for missing entries
Sigma_bar = isnan(Obs);
Sigma_bar = tensor(Sigma_bar,[nl nh nw]);

% nan to 0 and use it for Alg.
Obs(isnan(Obs)) = 0; 
D = tensor(Obs,[nl nh nw]);

Size = numel(double(Obs));


%% RTC
lambda = 1.4054;  % this parameter can be tuned

fprintf('lambda=%1.4f\n',lambda);

tic;
[Lhat,Shat, Ohat,iter] = inexact_alm_rmc21(D,Sigma_bar,lambda, 1e-7, 1000);
time = toc;

Spar = sum(sum(sum(double(Shat) >= 1e-5)))/Size;  % percentage of non-zero elements in S
disp(['Estimated Sparcity: ', num2str(Spar)])
% T = hosvd(Lhat,0.003);

%%%% examing results %%%
Shat_M = double(tenmat(Shat,1));
dates_ano = any(Shat_M >= 1e-3);   % test if there's any nonzero elements in a column
% sum(dates_ano)
dates_index = find(dates_ano);     % index of the column, the daylight shift is dealt with in python

% Lhat_M = double(tenmat(Lhat,1));
% dates_zero = all(Lhat_M <= 1e-5);
% sum(dates_ano)

dates_ano0 = [dates_ano(1:1659),dates_ano(1661:2856)];  % delete the additional entry for daylight time

% plot in columns
figure()
h = stem(dates_ano);
set(h, 'Marker', 'none');
% 
% % plot in blocks 
% Shat_3D = Shat.data;
% ST = any(Shat_3D >= 1e-3);
% ST = squeeze(ST);
% STI = ~ST;
% figure()
% imagesc([1,16],[1,168],STI)
% colormap('gray')
% 

% 
% Lhat_M = double(tenmat(Lhat,1));
% Lhat_M( :, ~any(Lhat_M,1) ) = [];  %delete columns with all zero
% L_mean = mean(Lhat_M,2);
% csvwrite('L_mean.txt',L_mean)
% mean & std

Lhat_3D = Lhat.data;
Lhat_3D(Lhat_3D < 1e-3) = NaN;
L_mean = nanmean(Lhat_3D,3); % infered mean by alg
L_std = nanstd(Lhat_3D,0,3);

L_mean0 = nanmean(Obs_Nan,3); % orig. mean
%% devi for mean;
I = 77 ;   % which hour to see;
I_hour = mod(I,168);  % which hour of the week
dev = (Obs_Nan(:,I) - L_mean(:,I_hour)) ./ L_std(:,I_hour);

figure()
plot(dev)
title( 'times of std from mean for all road segments')
%% plot
figure()
plot(double(D_plot(1,4,:)))
hold on

plot(double(D_plot(1,112,:)))
plot(double(D_plot(1,156,:)))
box
plot(double(D_plot(12,4,:)))
plot(double(D_plot(12,112,:)))
plot(double(D_plot(12,156,:)))

plot(double(D_plot(112,4,:)))
plot(double(D_plot(112,112,:)))
plot(double(D_plot(112,156,:)))

hold off


figure()
plot(double(D_plot(1,12,:)))
hold on
for i = 2:200
    
    plot(double(D_plot(i,12,:)))
end
hold off


figure()
plot(double(D_plot(1,:,4)))
hold on
for i = 2:200
    
    plot(double(D_plot(i,:,4)))
end
hold off

% %% curve for different lambda
% 
% n_lambda = 80;
% 
% % sweep
% Trank = zeros(n_lambda,1);  % Tucker rank an approximate
% Spar = zeros(n_lambda,1);  % Sparsity;
% err_S = zeros(n_lambda,1);  % Err of Sparse tensor vs. True S
% err_L = zeros(n_lambda,1);  % Err of Lowrank tensor vs. True L
% 
% Lambdas = linspace(0.1,120,n_lambda);      % all Lambdas tested
% 
% 
% 
% parfor j = 1:n_lambda
%     lambda = Lambdas(j);
%     Lambdas(j) = lambda;
%     [Lhat ,Shat,~] = inexact_alm_rpca21_v4(D, lambda, 1e-7, 1000);
%     try
%         T = hosvd(Lhat,0.3);
%         Trank(j) = sum(size(T.core));
%         Spar(j) = sum(sum(sum(double(Shat) >= 1e-5)))/Size;  % each sum operates along one array so call 3 times.
%     catch ex
%         fprintf('fail hosvd')
%     end
%     
%     fprintf('lambda=%1.4f\n',lambda);
%     fprintf('rank=%d;Spar = %1.4f\n',Trank(j),Spar(j));
% 
% end
% 
% %% plot sparsity & rank
% % series = step:step:(step*n_lambda);
% 
% figure();
% plot(Spar,Trank);
% xlabel('sparsity');
% ylabel('Trank');
% 
% 
% figure();
% plot(Lambdas,Trank );
% yyaxis right
% ylabel('Sparsity');
% hold on
% plot(Lambdas,Spar );
% xlabel('lambda');
% yyaxis left
% ylabel('Tucker Rank');
% 
% hold off
