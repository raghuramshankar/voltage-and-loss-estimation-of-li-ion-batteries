clear;
clc;
close all;

%% load variables
load loss_comparison.mat;

%% RMSE
rmse_loss = 1000*sqrt(mean((loss_phy(1:1221) - l4').^2))             % [mW]

%% plot
figure(1);
plot(l4, 'r', 'LineWidth', 2);
hold on;
plot(loss_phy(1:1221), 'g', 'LineWidth', 2);
hold on;
plot(loss_rc(1:1221), 'b--', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Loss [W]');
title('Loss Comparison - MATLAB physical, COMSOL physical, RC');
legend('MATLAB', 'COMSOL', 'RC');