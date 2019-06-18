% this code loads simulation data and draws grayscale uimage of possibility
% of success for given ibservationratio and rank

addpath PROPACK/tensor_toolbox-master
addpath ../SimulationData

clear;
load sim_block_v3.mat    % data simulated from draw_block.m

x1 =  size(Precision,2); % ratio
x2 = size(Precision,3); % rank

Pre = Precision>= 0.999;
Rec = Recall>= 0.999;

Suc = Pre & Rec;   % define success as both precisin and recall larger than 0.99
SucM = mean(Suc);  % calculate the posibility of success across all trials.
SucM = reshape(SucM, x1,x2);

SucM_s = SucM(:,1:19);

% plot
figure()
imagesc([2,20],[0.3,1],SucM_s)
colormap('gray')
bar = colorbar;
bar.Label.String = 'Rate of success';
xlabel('Rank c');
ylabel('Observation ratio');
set(gca, 'fontsize', 18)
title('Rate of Successful Outlier Identification')

% save
print('../Figure/block','-depsc')



