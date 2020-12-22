clear;
clc;
close all;

%% load data
load('ocv.mat');
ocv = ocv';
soc = (linspace(0.9, 0.1, size(ocv, 1)))';

%% postprocessing
figure(1);
plot(soc, ocv, 'r', 'LineWidth', 2);
xlabel('SOC');
ylabel('OCV [V]');
title('OCV Curve');
legend('OCV', 'Location', 'South West');