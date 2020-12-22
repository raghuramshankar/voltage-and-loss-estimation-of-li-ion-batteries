close all
clear
clc

%% initial condition
% 4.16165 V -> 90% SOC
% 3.33541 V -> 10% SOC

% %% pulse 0
% soc_init_neg = 0.9; 
% soc_init_pos = 0.1; 

% % pulse 1
% soc_init_neg = 0.8155;  
% soc_init_pos = 0.15; 

% % pulse 1
% soc_init_neg = 0.82;  
% soc_init_pos = 0.14; 

% % pulse 2
% soc_init_neg = 0.73056; 
% soc_init_pos = 0.2024;

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
% % pulse 6
% soc_init_neg = 0.392; 
% soc_init_pos = 0.407;

% %% pulse 7
% soc_init_neg = 0.30757; 
% soc_init_pos = 0.458;

% %% pulse 8
% soc_init_neg = 0.223; 
% soc_init_pos = 0.509;
% 
%% pulse 9
soc_init_neg = 0.13836; 
soc_init_pos = 0.56;

% %% pulse 10
% soc_init_neg = 0.05375; 
% soc_init_pos = 0.6114; 

%% time and voltage limations
time_max = 3600; % seconds
% time_max = time_max - 1;
time_step = 1;
V_max = 4.2;
V_min = 0.0;

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
% clc;
% I = 1e-5*ones(floor(4*s/10), 1);
% I = [I; 0.25*c_1];
% I = [I; 0.5*c_1];
% I = [I; 0.75*c_1];
% I = [I; 1*c_1*ones(floor(s/10) - 6, 1)];
% c_1 = 1e-5;
% I = [I; 0.75*c_1];
% I = [I; 0.5*c_1];
% I = [I; 0.25*c_1];
% I = [I; 1e-5*ones(floor(5*s/10), 1)];

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
I = 1e-5;

Setup_parameters;
SpatialDiscretization;
Initialize;

% initialize data storage
time = 0:time_step:time_max; % only valid for a constant time step
Ah_time = zeros(time_max/time_step+1,1);
V_batt_time = zeros(time_max/time_step+1,1);V_batt_time(1) = phis_eq_pos(end)-phis_eq_neg(1);
c_s_surface_to_compare = zeros(time_max/time_step+1,n_mesh_neg+n_mesh_pos+2);c_s_surface_to_compare(1,:) = [c_s_surface_neg c_s_surface_pos];
c_l_to_compare=zeros(time_max/time_step+1,n_mesh_neg+n_mesh_sep+n_mesh_pos+3);c_l_to_compare(1,:) = c_l;
time_now = 0;
k=1;
tic

% I = 1e-5;
% time_max = 3600*1/10; % seconds
% time_max = 1;
% Simulation_loop;

I = 16.4;
time_max = 3600*1/10; % seconds
CurrentStep_prep;
Simulation_loop;

I = 1e-5;
time_max = 3600*10/10; % seconds
CurrentStep_prep
Simulation_loop;
toc

%% postprocessing
Area = Area * 1e4; % convert m^2 to cm^2
Ah_time = Ah_time ./ 1e3; % convert mAh to Ah
Ah_time = Ah_time .* Area; % convert mAh/cm^2 to Ah

I = [1e-5; 16.4*ones(3600*1/10 + 1, 1); 1e-5*ones(3600*9/10 - 2, 1)]; 

% figure(1)
% plot(Ah_time,V_batt_time,'LineWidth',2);
% 
% hold on
% plot(Ah_time,V_batt_time,'LineWidth',2);
% grid on
% xlabel('Charge (Ah)')
% ylabel('Voltage (V)')

figure(2)
plot(time, V_batt_time, 'LineWidth', 2);
grid on
xlabel('Time [s]')
ylabel('Voltage (V)')