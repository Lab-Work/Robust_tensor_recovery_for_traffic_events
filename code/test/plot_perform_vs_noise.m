%this code plots curve of results v.s observation ratio from tensor
%completion with varying corruption ratio and the same rank (5,5,5)

addpath ../SimulationData

% load first data set
load sim_rmc_0.05_70_5.mat
figure(3)
plot(Ratio_D, err_L_avg,'LineWidth', 2 );
xlabel('observation ratio');
ylabel('RE');
title({'RE v.s. observation ratio','for different corruption ratios'});
ylim([0,1]);
set(gca, 'fontsize', 18)
hold on;

figure(4);
plot(Ratio_D, Recall_avg, 'LineWidth', 2 );
xlabel('observation ratio');
ylabel('Recall');
title({'Recall v.s. observation ratio','for different corruption ratios'});
ylim([0,1.1]);
set(gca, 'fontsize', 18)
hold on

figure(5);
plot(Ratio_D, Precision_avg, 'LineWidth', 2 );
xlabel('observation ratio');
ylabel('Precision');
title({'Precision v.s. observation ratio','for different corruption ratios'});
ylim([0,1.1]);
set(gca, 'fontsize', 18)
hold on


figure(6);
plot(Ratio_D, Iter_avg, 'LineWidth', 2 );
xlabel('observation ratio');
ylabel('# Iteration');
title({'# Iteration v.s. observation ratio','for different corruption ratios'});
set(gca, 'fontsize', 18)
ylim([0,42]);
hold on

% load second data set
load sim_rmc_0.1_70_5(1).mat
figure(3)
plot(Ratio_D, err_L_avg,'LineWidth', 2 );

figure(4);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );

figure(5);
plot(Ratio_D, Precision_avg,'LineWidth', 2 );

figure(6);
plot(Ratio_D, Iter_avg,'LineWidth', 2 );


% load last set
load sim_rmc_0.2_70_5(1).mat
figure(3)
plot(Ratio_D, err_L_avg,'LineWidth', 2 );
lgd = legend('\gamma = 0.05','\gamma = 0.1','\gamma = 0.2');
title(lgd, 'corruption ratio')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/err_L_cor','-depsc')

figure(4);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );
lgd = legend('\gamma = 0.05','\gamma = 0.1','\gamma = 0.2');
title(lgd, 'corruption ratio')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/recall_cor','-depsc')


figure(5);
plot(Ratio_D, Precision_avg,'LineWidth', 2 );
lgd = legend('\gamma = 0.05','\gamma = 0.1','\gamma = 0.2');
title(lgd, 'corruption ratio')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/precison_cor','-depsc')


figure(6);
plot(Ratio_D, Iter_avg,'LineWidth', 2 );
lgd = legend('\gamma = 0.05','\gamma = 0.1','\gamma = 0.2');
title(lgd, 'corruption ratio')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/iter_cor','-depsc')
