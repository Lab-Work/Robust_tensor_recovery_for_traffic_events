% this code compares the performance of rpca using l1 and l2,1 norm.

addpath PROPACK/tensor_toolbox-master
addpath ../SimulationData

% load first data set
load sim_rpca3D(1).mat
figure(3)
plot(Ratio_s, err_L_avg,'LineWidth', 2 );
xlabel('corruption ratio');
ylabel('RE');
title({'RE v.s. corruption ratio','for l_1 and l_{2,1} norm regularizations'});
% ylim([0,1]);
set(gca, 'fontsize', 18)
hold on;

figure(4);
plot(Ratio_s, Recall_avg, 'LineWidth', 2 );
xlabel('corruption ratio');
ylabel('Recall');
title({'Recall v.s. corruption ratio','for l_1 and l_{2,1} norm regularizations'});
ylim([0,1.1]);
set(gca, 'fontsize', 18)
hold on

figure(5);
plot(Ratio_s, Precision_avg, 'LineWidth', 2 );
xlabel('corruption ratio');
ylabel('Precision');
title({'Precision v.s. corruption ratio','for l_1 and l_{2,1} norm regularizations'});
ylim([0,1.1]);
set(gca, 'fontsize', 18)
hold on


figure(6);
plot(Ratio_s, Iter_avg, 'LineWidth', 2 );
xlabel('corruption ratio');
ylabel('# Iteration');
title({'# Iteration v.s. corruption ratio','for l_1 and l_{2,1} norm regularizations'});
set(gca, 'fontsize', 18)
ylim([0,42]);
hold on



% load second set
load sim_rpca.mat
figure(3)
plot(Ratio_s, err_L_avg,'LineWidth', 2 );
lgd = legend('l_1 norm','l_{2,1} norm');
title(lgd, 'regularization')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/err_noise_compare','-depsc')

figure(4);
plot(Ratio_s, Recall_avg,'LineWidth', 2 );
lgd = legend('l_1 norm','l_{2,1} norm');
title(lgd, 'regularization')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/recall_noise_compare','-depsc')


figure(5);
plot(Ratio_s, Precision_avg,'LineWidth', 2 );
lgd = legend('l_1 norm','l_{2,1} norm');
title(lgd, 'regularization')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/precison_noise_compare','-depsc')


figure(6);
plot(Ratio_s, Iter_avg,'LineWidth', 2 );
lgd = legend('l_1 norm','l_{2,1} norm');
title(lgd, 'regularization')
set(lgd, 'Location', 'Best')
hold off
print('../Figure/iter_noise_compare','-depsc')
