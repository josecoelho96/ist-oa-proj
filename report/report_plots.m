close all;
clear;

load workspace_10-12-17_23-21-33.mat
clearvars -except results_mse
results_mse_part1 = results_mse;

load workspace_13-12-17_11-45-39.mat
clearvars -except results_mse_part1 results_mse
results_mse_part2 = results_mse;
results_mse = [ results_mse_part1  results_mse_part2'];

SNR = [5 10 15 20 25]; % SNR wanted
semilogy(results_mse, '.-', 'MarkerSize',20, 'LineWidth', 1.5)
title('MSE variation with SNR')
xlabel('SNR [dB]')
ylabel('MSE')
legend('LS (GA)', 'LS', 'L_1', 'P_1', 'P_2(1)', 'MS(1)', 'Location', 'southwest');
ax = gca;
ax.XTick = [1 2 3 4 5];
ax.XTickLabel = SNR;
grid on;