% Example script to simulate a discharge process
close all
clear
clc

% initial condition
soc_init_neg = 0.9;
soc_init_pos = 0.1;

% time and voltage limations
time_max = 1e6; % seconds
V_max = 4.2; % cut-over voltage
V_min = 3; % cut-off voltage

time_step = 10;
order = 4;
I = 20; %A/m^2

Setup_parameters_b; % The cell geometry, material properties, meshing are defined here
SpatialDiscretization_b; % create the spatial discretization
Initialize_b;

tic
Simulation_loop_b;
toc

%% plot
figure
hold on
plot(Ah_time,V_batt_time,'LineWidth',2);
grid on
xlabel('Charge (mAh/cm^2)')
ylabel('Voltage (V)')

