%this code plots 3 curve of results v.s. observation ratio from tensor
%completion with varying tensor rank c and the same tensor corruption ratio
%at 0.1

addpath PROPACK/tensor_toolbox-master
addpath ../../SimulationData

% load first data set
load sim_rmc_0.1_70_8(1).mat
figure(1)
plot(Ratio_D, err_L_avg,'LineWidth', 2  );
xlabel('observation ratio');
ylabel('RES');
title({'RES v.s. observation ratio','for different tucker ranks'});
ylim([0,1.1]);
set(gca, 'fontsize', 18)
hold on;

figure(2);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );
xlabel('observation ratio');
ylabel('Recall rate');
title({'Recall v.s. observation ratio','for different tucker ranks'});
ylim([0,1]);
set(gca, 'fontsize', 18)
hold on

figure(3);
plot(Ratio_D, Precision_avg, 'LineWidth', 2 );
xlabel('observation ratio');
ylabel('Precision');
title({'Precision v.s. observation ratio','for different tucker ranks'});
ylim([0,1.1]);
set(gca, 'fontsize', 18)
hold on


figure(4);
plot(Ratio_D, Iter_avg, 'LineWidth', 2 );
xlabel('observation ratio');
ylabel('# Iteration');
title({'# Iteration v.s. observation ratio','for different tucker ranks'});
set(gca, 'fontsize', 18)
ylim([1,45]);
hold on

% load second data set
load sim_rmc_0.1_70_5(1).mat
figure(1)
plot(Ratio_D, err_L_avg,'LineWidth', 2 );

figure(2);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );

figure(3);
plot(Ratio_D, Precision_avg,'LineWidth', 2 );

figure(4);
plot(Ratio_D, Iter_avg,'LineWidth', 2 );


% load last data set
load sim_rmc_0.1_70_2(1).mat
figure(1)
plot(Ratio_D, err_L_avg,'LineWidth', 2  );
lgd = legend('(8,8,8)','(5,5,5)','(2,2,2)');
title(lgd, 'tucker rank')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/err_rank_3','-depsc')

figure(2);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );
lgd = legend('(8,8,8)','(5,5,5)','(2,2,2)');
title(lgd, 'tucker rank')
set(lgd, 'Location', 'Best')
ylim([0,1.1]);
hold off
print('../Figure/recall_rank_3','-depsc')


figure(3);
plot(Ratio_D, Precision_avg,'LineWidth', 2 );
lgd = legend('(8,8,8)','(5,5,5)','(2,2,2)');
title(lgd, 'tucker rank')
set(lgd, 'Location', 'Best')
ylim([0,1.1]);
hold off
print('../Figure/precison_rank_3','-depsc')


figure(4);
plot(Ratio_D, Iter_avg,'LineWidth', 2 );
lgd = legend('(8,8,8)','(5,5,5)','(2,2,2)');
title(lgd, 'tucker rank')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/iter_rank_3','-depsc')
