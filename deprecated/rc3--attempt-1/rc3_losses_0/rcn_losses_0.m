clear;
clc;
close all;

%% load data
load('rc_zero.mat');
load('rc3_0.mat');

%% variables
s = length(time);
t_start = floor(4*s/10);
t_end = floor(5*s/10);
t_pulse = t_end - t_start;
I = max(I);

%% calculate losses
area_zero = trapz(time(t_start:t_end), v_t_zero(t_start:t_end));
area_rc = area_zero - trapz(time(t_start:t_end), v_t(t_start:t_end));
area_phy = area_zero - trapz(time(t_start:t_end), V_batt_time(t_start:t_end));

loss_rc = I * area_rc / t_pulse                                 % [W]
loss_phy = I * area_phy / t_pulse                               % [W]   

r_rc = loss_rc / (I^2)                                          % [ohm] 
r_phy = loss_phy / (I^2)                                        % [ohm]

%% postprcoessing
figure(1);
plot(time, v_t, 'g', 'LineWidth', 2);
hold on;
plot(time, v_t_zero, 'r', 'LineWidth', 2);
hold on;
plot(time, V_batt_time, 'b', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Terminal Voltage [V]');
title('Loss Area');
legend('RC3', 'No RC', 'Physical');
