%this code plots 5 curve of results v.s. observation ratio from tensor
%completion with varying tensor rank c and the same tensor corruption ratio
%at 0.1

addpath PROPACK/tensor_toolbox-master
addpath ../SimulationData

% load first data set
load sim_rmc_0.1_70_10.mat
figure(1)
plot(Ratio_D, err_L_avg );
xlabel('observstion ratio');
ylabel('relative error of L');
title({'Relative error v.s. observation ratio','for different tucker rank'});
ylim([0,1]);
set(gca, 'fontsize', 18)
hold on;

figure(2);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );
xlabel('observstion ratio');
ylabel('Recall rate');
title({'Recall v.s. observation ratio','for different tucker rank'});
ylim([0,1]);
set(gca, 'fontsize', 18)
hold on

figure(3);
plot(Ratio_D, Precision_avg, 'LineWidth', 2 );
xlabel('observstion ratio');
ylabel('Precision');
title({'Precision v.s. observation ratio','for different tucker rank'});
ylim([0,1]);
set(gca, 'fontsize', 18)
hold on


figure(4);
plot(Ratio_D, Iter_avg, 'LineWidth', 2 );
xlabel('observstion ratio');
ylabel('# Iteration');
title({'# Iteration v.s. observation ratio','for different tucker rank'});
set(gca, 'fontsize', 18)
ylim([13,37]);
hold on

% load second data set
load sim_rmc_0.1_70_8.mat
figure(1)
plot(Ratio_D, err_L_avg,'LineWidth', 2 );

figure(2);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );

figure(3);
plot(Ratio_D, Precision_avg,'LineWidth', 2 );

figure(4);
plot(Ratio_D, Iter_avg,'LineWidth', 2 );

% % load data set
% load sim_rmc_0.1_70_5.mat
% figure(1)
% plot(Ratio_D, err_L_avg,'LineWidth', 2 );
% 
% figure(2);
% plot(Ratio_D, Recall_avg,'LineWidth', 2 );
% 
% figure(3);
% plot(Ratio_D, Precision_avg,'LineWidth', 2 );
% 
% figure(4);
% plot(Ratio_D, Iter_avg,'LineWidth', 2 );


% load data set
load sim_rmc_0.1_70_5.mat
figure(1)
plot(Ratio_D, err_L_avg,'LineWidth', 2 );

figure(2);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );

figure(3);
plot(Ratio_D, Precision_avg,'LineWidth', 2 );

figure(4);
plot(Ratio_D, Iter_avg,'LineWidth', 2 );


% load last data set
load sim_rmc_0.1_70_2.mat
figure(1)
plot(Ratio_D, err_L_avg );
lgd = legend('c=10','c=8','c=5','c=2');
title(lgd, 'tucker rank is (c,c,c)')
set(lgd, 'Location', 'Best')
hold off
% print('../Figure/err_rank_5','-depsc')

figure(2);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );
lgd = legend('c=10','c=8','c=5','c=2');
title(lgd, 'tucker rank is (c,c,c)')
set(lgd, 'Location', 'Best')
hold off
% print('../Figure/recall_rank_5','-depsc')


figure(3);
plot(Ratio_D, Precision_avg,'LineWidth', 2 );
lgd = legend('c=10','c=8','c=5','c=2');
title(lgd, 'corruption ratio')
set(lgd, 'Location', 'Best')
hold off
% print('../Figure/precison_rank_5','-depsc')


figure(4);
plot(Ratio_D, Iter_avg,'LineWidth', 2 );
% lgd = legend('0.1','0.3','0.3');
lgd = legend('c=10','c=8','c=5','c=2');
title(lgd, 'corruption ratio')
set(lgd, 'Location', 'Best')
hold off
% print('../Figure/iter_rank_5','-depsc')
