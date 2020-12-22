close all
clear
clc

%% initial condition
% 4.16165 V -> 90% SOC
% 3.33541 V -> 10% SOC

%% pulse 0
soc_init_neg = 0.9; 
soc_init_pos = 0.1; 

% % pulse 1
% soc_init_neg = 0.8155;  
% soc_init_pos = 0.15; 

% %% pulse 2
% soc_init_neg = 0.73056; 
% soc_init_pos = 0.2024;
% 
% %% pulse 3
% soc_init_neg = 0.646; 
% soc_init_pos = 0.2535; 
% 
% %% pulse 4
% soc_init_neg = 0.5615; 
% soc_init_pos = 0.3045;

% %% pulse 5
% soc_init_neg = 0.4765; 
% soc_init_pos = 0.3558;
% 
% %% pulse 6
% soc_init_neg = 0.392; 
% soc_init_pos = 0.407;

% %% pulse 7
% soc_init_neg = 0.30757; 
% soc_init_pos = 0.458;
% 
% %% pulse 8
% soc_init_neg = 0.223; 
% soc_init_pos = 0.509;
% 
% %% pulse 9
% soc_init_neg = 0.13836; 
% soc_init_pos = 0.56;
% 
% %% pulse 10
% soc_init_neg = 0.05375; 
% soc_init_pos = 0.6114; 

%% time and voltage limations
time_max = 3600; % seconds
% time_max = time_max - 1;
time_step = 1;
V_max = 4.2;
V_min = 3.0;

order = 4;
time = 0:time_step:time_max;
s = length(time);
soc = 0.0:1/s:1.0;
c_1 = 16.4; % A/m^2

%% constant discharge
% clc;
% I = c_1/20*ones(s,1); %A/m^2
% I = -c_1/20*ones(s,1); %A/m^2

%% pulse discharge
% % clc;
% I = 1e-5*ones(floor(4*s/10), 1);
% I = [I; 0.25*c_1];
% I = [I; 0.5*c_1];
% I = [I; 0.75*c_1];
% I = [I; 1*c_1*ones(floor(s/10) - 6, 1)];
% I = [I; 0.75*c_1];
% I = [I; 0.5*c_1];
% I = [I; 0.25*c_1];
% I = [I; 1e-5*ones(floor(5*s/10), 1)];

%% nedc driving cycle
clc;
load('nedc.mat');
I_0 = out.current_NEDC.Data(:, 1);
I = 1e-5 + I_0;
time = 0:time_step:length(I);

%% run scripts
Setup_parameters;
SpatialDiscretization;
Initialize;

tic
Simulation_loop;
toc

%% postprocessing
Area = Area * 1e4; % convert m^2 to cm^2
Ah_time = Ah_time ./ 1e3; % convert mAh to Ah
Ah_time = Ah_time .* Area; % convert mAh/cm^2 to Ah

%% plot
figure(1)
plot(time, I,'LineWidth',2);
grid on
xlabel('Time (s)')
ylabel('Curent (A)')

figure(2)
plot(time, V_batt_time, 'LineWidth', 2);
grid on
xlabel('Time (s)')
ylabel('Voltage (V)')