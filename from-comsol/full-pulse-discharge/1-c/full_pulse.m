clear;
clc;
close all;

%% load data
load('current.txt')
load('voltage.txt')
I = current(:, 2);
V_batt_time = voltage(:, 2);
total_time = voltage(:, 1);
time_step = 1;
% soc(interp1(ocv, 1:length(ocv), V_batt_time(1), 'nearest'))

%% pulse 0
V_batt_time = voltage(1:3604, 2);
I = current(1:3604, 2);
time = voltage(1:3604, 1);
soc_init_neg = 0.9;
soc_init_pos = 0.1;

%% pulse 1
V_batt_time = voltage(3604:7206, 2);
I = current(3604:7206, 2);
time = voltage(3604:7206, 1);
soc_init_neg = 0.8155;
soc_init_pos = 0.15;

%% pulse 2
V_batt_time = voltage(7206:10806, 2);
I = current(7206:10806, 2);
time = voltage(7206:10806, 1);
soc_init_neg = 0.73056;
soc_init_pos = 0.2024;

%% pulse 3
V_batt_time = voltage(10806:14406, 2);
I = current(10806:14406, 2);
time = voltage(10806:14406, 1);
soc_init_neg = 0.6460;
soc_init_pos = 0.2535;

%% pulse 4
V_batt_time = voltage(14406:18005, 2);
I = current(14406:18005, 2);
time = voltage(14406:18005, 1);
soc_init_neg = 0.5615;
soc_init_pos = 0.3045;

%% pulse 5
V_batt_time = voltage(18005:21606, 2);
I = current(18005:21606, 2);
time = voltage(18005:21606, 1);
soc_init_neg = 0.4765;
soc_init_pos = 0.3558;

%% pulse 6
V_batt_time = voltage(21606:25206, 2);
I = current(21606:25206, 2);
time = voltage(21606:25206, 1);
soc_init_neg = 0.3920;
soc_init_pos = 0.4070;

%% pulse 7
V_batt_time = voltage(25206:28805, 2);
I = current(25206:28805, 2);
time = voltage(25206:28805, 1);
soc_init_neg = 0.3076;
soc_init_pos = 0.4580;

%% pulse 8
V_batt_time = voltage(28805:32402, 2);
I = current(28805:32402, 2);
time = voltage(28805:32402, 1);
soc_init_neg = 0.2230;
soc_init_pos = 0.5090;

%% pulse 9
V_batt_time = voltage(32402:36002, 2);
I = current(32402:36002, 2);
time = voltage(32402:36002, 1);
soc_init_neg = 0.1384;
soc_init_pos = 0.5600;