%this code compares the performance of robust tensor completion using l1 and l2,1 norm.
addpath PROPACK/tensor_toolbox-master
addpath ../SimulationData

% load first data set
load sim_rmc3D_0.1_70_5.mat
figure(3)
plot(Ratio_D, err_L_avg,'LineWidth', 2 );
xlabel('observation ratio');
ylabel('RE');
title({'RE v.s. observation ratio','for l_1 and l_{2,1} norm regularizations'});
ylim([0,1]);
set(gca, 'fontsize', 18)
hold on;

figure(4);
plot(Ratio_D, Recall_avg, 'LineWidth', 2 );
xlabel('observation ratio');
ylabel('Recall');
title({'Recall v.s. observation ratio','for l_1 and l_{2,1} norm regularizations'});
ylim([0,1.1]);
set(gca, 'fontsize', 18)
hold on

figure(5);
plot(Ratio_D, Precision_avg, 'LineWidth', 2 );
xlabel('observation ratio');
ylabel('Precision');
title({'Precision v.s. observation ratio','for l_1 and l_{2,1} norm regularizations'});
ylim([0,1.1])
set(gca, 'fontsize', 18)
hold on


figure(6);
plot(Ratio_D, Iter_avg, 'LineWidth', 2 );
xlabel('observation ratio');
ylabel('# Iterations');
title({'# Iteration v.s. observation ratio','for l_1 and l_{2,1} norm regularizations'});
set(gca, 'fontsize', 18)
ylim([0,42]);
hold on


% load last set
load sim_rmc_0.1_70_5(1).mat
figure(3)
plot(Ratio_D, err_L_avg,'LineWidth', 2 );
lgd = legend('l_1 norm','l_{2,1} norm');
title(lgd, 'regularization')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/err_L_cor_rmccompare','-depsc')

figure(4);
plot(Ratio_D, Recall_avg,'LineWidth', 2 );
lgd = legend('l_1 norm','l_{2,1} norm');
title(lgd, 'regularization')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/recall_cor_rmccompare','-depsc')


figure(5);
plot(Ratio_D, Precision_avg,'LineWidth', 2 );
lgd = legend('l_1 norm','l_{2,1} norm');
title(lgd, 'regularization')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/precison_cor_rmccompare','-depsc')


figure(6);
plot(Ratio_D, Iter_avg,'LineWidth', 2 );
lgd = legend('l_1 norm','l_{2,1} norm');
title(lgd, 'regularization')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/iter_cor_rmccompare','-depsc')
