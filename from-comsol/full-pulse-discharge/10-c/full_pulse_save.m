clear;
clc;
% close all;

%% load data
load PulseDischargeData_10C.mat
load curr1.mat

%% pulse 0
V_batt_time = Voltage_10C(1:3603);
I = curr1(1:3603, 2);
time = time(1:3603);

%% pulse 1
V_batt_time = Voltage_10C(3603:7205);
I = curr1(3603:7205, 2);
time = time(3603:7205);

%% pulse 2
V_batt_time = Voltage_10C(7205:10805);
I = curr1(7205:10805, 2);
time = time(7205:10805);

%% pulse 3
V_batt_time = Voltage_10C(10805:14403);
I = curr1(10805:14403, 2);
time = time(10805:14403);

%% pulse 4
V_batt_time = Voltage_10C(14403:18002);
I = curr1(14403:18002, 2);
time = time(14403:18002);

%% pulse 5
V_batt_time = Voltage_10C(18002:21601);
I = curr1(18002:21601, 2);
time = time(18002:21601);

%% pulse 6
V_batt_time = Voltage_10C(21601:25200);
I = curr1(21601:25200, 2);
time = time(21601:25200);

%% pulse 7
V_batt_time = Voltage_10C(25200:28805);
I = curr1(25200:28805, 2);
time = time(25200:28805);

%% load data
load PulseDischargeData_10C_lastfew.mat
load curr2.mat
time = time(1:10817);

%% pulse 8
V_batt_time = voltage_10c(1:3604);
time = time(1:3604);
time = time + 28805;
I = curr2(1:3604, 2);

%% pulse 9
V_batt_time = voltage_10c(3604:7205);
time = time(3604:7205);
time = time + 28805;
I = curr2(3604:7205, 2);

%% pulse 10
V_batt_time = voltage_10c(7205:10805);
time = time(7205:10805);
time = time + 28805;
I = curr2(7205:10805, 2);
