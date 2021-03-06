clear;
clc;
close all;

%% load data
load('current.txt')
load('voltage.txt')
I = current(:, 2);
V_batt_time = voltage(:, 2);
time = voltage(:, 1);
time_step = 1;
% soc(interp1(ocv, 1:length(ocv), V_batt_time(1), 'nearest'))