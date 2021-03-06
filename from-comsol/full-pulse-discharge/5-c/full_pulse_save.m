clear;
clc;
% close all;

%% load data
load PulseDischargeData_5C.mat
load curr1.mat

%% pulse 0
V_batt_time = Voltage_5C(1:3605);
I = curr1(1:3605, 2);
time = time(1:3605);

%% pulse 1
V_batt_time = Voltage_5C(3605:7205);
I = curr1(3605:7205, 2);
time = time(3605:7205);

%% pulse 2
V_batt_time = Voltage_5C(7205:10805);
I = curr1(7205:10805, 2);
time = time(7205:10805);

%% pulse 3
V_batt_time = Voltage_5C(10805:14404);
I = curr1(10805:14404, 2);
time = time(10805:14404);

%% pulse 4
V_batt_time = Voltage_5C(14404:18003);
I = curr1(14404:18003, 2);
time = time(14404:18003);

%% pulse 5
V_batt_time = Voltage_5C(18003:21602);
I = curr1(18003:21602, 2);
time = time(18003:21602);

%% pulse 6
V_batt_time = Voltage_5C(21602:25205);
I = curr1(21602:25205, 2);
time = time(21602:25205);

%% pulse 7
V_batt_time = Voltage_5C(25205:28800);
I = curr1(25205:28800, 2);
time = time(25205:28800);

%% pulse 8
V_batt_time = Voltage_5C(28800:32398);
time = time(28800:32398);
I = curr1(28800:32398, 2);
