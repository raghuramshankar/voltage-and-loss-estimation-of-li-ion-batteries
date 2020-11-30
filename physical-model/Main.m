close all
clear
clc

%% initial condition
soc_init_neg = 0.9;
soc_init_pos = 0.1;

%% time and voltage limations
time_max = 1000; % seconds
time_step = 1;
V_max = 4.2; % cut-over voltage
V_min = 3; % cut-off voltage

order = 4;
time = 0:time_step:time_max;
s = size(time, 2);
soc = 0.1:1/s:0.9;

%% constant discharge
% clc;
I = 20*ones(s,1); %A/m^2

% %% pulse discharge
% clc;
I = 1e-5*ones(floor(4*s/10), 1);
I = [I; 10*ones(floor(s/10), 1)];
I = [I; 1e-5*ones(floor(5*s/10), 1)];

%% nedc driving cycle
% clc;
% load('nedc.mat');
% I_0 = out.current_NEDC.Data(:, 1);
% for i = 1:size(time, 2)
%     if i > size(I_0, 1)
%         i = i - size(I_0, 1);
%     end
%     I(i) = I_0(i) + 1e-5;
% end

%% run scripts
clc
Setup_parameters; % The cell geometry, material properties, meshing are defined here
SpatialDiscretization; % create the spatial discretization
Initialize;

tic
Simulation_loop;
toc

%% plot
figure(1)
plot(Ah_time,V_batt_time,'LineWidth',2);
grid on
xlabel('Charge (mAh/cm^2)')
ylabel('Voltage (V)')

figure(2)
plot(time, V_batt_time, 'LineWidth', 2);
grid on
xlabel('Time [s]')
ylabel('Voltage (V)')