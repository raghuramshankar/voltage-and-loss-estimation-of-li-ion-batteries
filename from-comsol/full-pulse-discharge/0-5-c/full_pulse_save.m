clear;
clc;
% close all;

%% load data
load PulseDischargeData_0.5C.mat
load curr1.mat

%% pulse 0
V_batt_time = Voltage(1:3603);
I = curr1(1:3603, 2);
time = time(1:3603);

%% pulse 1
V_batt_time = Voltage(3603:7205);
I = curr1(3603:7205, 2);
time = time(3603:7205);

%% pulse 2
V_batt_time = Voltage(7205:10805);
I = curr1(7205:10805, 2);
time = time(7205:10805);

%% pulse 3
V_batt_time = Voltage(10805:14405);
I = curr1(10805:14405, 2);
time = time(10805:14405);

%% pulse 4
V_batt_time = Voltage(14405:18005);
I = curr1(14405:18005, 2);
time = time(14405:18005);

%% pulse 5
V_batt_time = Voltage(18005:21605);
I = curr1(18005:21605, 2);
time = time(18005:21605);

%% pulse 6
V_batt_time = Voltage(21605:25205);
I = curr1(21605:25205, 2);
time = time(21605:25205);

%% pulse 7
V_batt_time = Voltage(25205:28805);
I = curr1(25205:28805, 2);
time = time(25205:28805);

%% pulse 8
V_batt_time = Voltage(28805:32404);
time = time(28805:32404);
I = curr1(28805:32404, 2);

%% pulse 9
V_batt_time = Voltage(32404:36005);
time = time(32404:36005);
I = curr1(32404:36005, 2);

%% pulse 10
V_batt_time = Voltage(36005:39606);
time = time(36005:39606);
I = curr1(36005:39606, 2);

%% pulse 11
V_batt_time = Voltage(39606:43206);
time = time(39606:43206);
I = curr1(39606:43206, 2);

%% pulse 12
V_batt_time = Voltage(43206:46806);
time = time(43206:46806);
I = curr1(43206:46806, 2);

%% pulse 13
V_batt_time = Voltage(46806:50406);
time = time(46806:50406);
I = curr1(46806:50406, 2);

%% pulse 14
V_batt_time = Voltage(50406:54006);
time = time(50406:54006);
I = curr1(50406:54006, 2);

%% pulse 15
V_batt_time = Voltage(54006:57605);
time = time(54006:57605);
I = curr1(54006:57605, 2);

%% pulse 16
V_batt_time = Voltage(57605:61193);
time = time(57605:61193);
I = curr1(57605:61193, 2);

%% pulse 17
V_batt_time = Voltage(61193:64802);
time = time(61193:64802);
I = curr1(61193:64802, 2);

%% pulse 18
V_batt_time = Voltage(64802:68397);
time = time(64802:68397);
I = curr1(64802:68397, 2);

%% pulse 19
V_batt_time = Voltage(68397:71996);
time = time(68397:71996);
I = curr1(68397:71996, 2);