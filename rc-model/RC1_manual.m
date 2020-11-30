clear;
clc;
close all;

%% load data
load('volt_vs_soc_1a.mat');
ocv = V_batt_time;
load('pulse_10a.mat')

%% time and voltage limitations
v_max = 4.2;                                % cut-over voltage
v_min = 3;                                  % cut-off voltage
s = size(time, 2);

%% variables
dt = 1;                                     % sample time [s]
A = 1*1e4;                                  % area of cell [cm^2]
% Q = Ah_time(end)*A*1e-3;                    % capacity [Ah]
Q = 1.6;
eta = 1;                                    % coulombic efficiency

i_r = zeros(s, 1);
v_c = zeros(s, 1);
% v_t = zeros(s, 1);
ind = zeros(s, 1);
v_t(1) = ocv(1);
soc = linspace(0.9, 0.1, s);
flag = 1;

%% manual fit
r0 = 9*1e-3;                               % [ohms]
r = 4*1e-3;                                % [ohms]
c = 70*1e3;                                % [farad]
% r = 0;
% c = 0;

%% pulse discharge
I = 1e-5*ones(floor(4*s/10), 1);
I = [I; 10*ones(floor(s/10), 1)];
I = [I; 1e-5*ones(floor(5*s/10), 1)];

%% compute all
f = exp(-dt/(r*c));
z0 = soc(1);
for k = 1:s-1
    if flag == 1
        z(k) = z0 - dt/(Q*3600)*eta*sum(I(1:k));
        i_r(k+1) = diag(f)*i_r(k) + (1-f)*I(k);
        v_c(k+1) = r.*i_r(k+1);
        ind(k) = interp1(soc, 1:length(soc), z(k), 'nearest');
        v_t(k+1) = ocv(ind(k)) - v_c(k) - I(k).*r0;
        if v_t(k) < v_min
            fprintf('Reached lower voltage limit.');
            flag = 0;
        end
        if v_t(k) > v_max
            fprintf('Reached upper voltage limit.');
            flag = 0;
        end 
    end
end

%% compute rms error
rmse = 1000*sqrt(mean((V_batt_time' - v_t).^2))     % [mV]

%% plot
figure(1);
plot(time, v_t, 'r', 'LineWidth', 2)
hold on;
plot(time, V_batt_time, 'b', 'LineWidth', 2)
xlabel('Time [s]')
ylabel('Terminal Voltage [V]')
title('RC1')
legend('RC1', 'Physical', 'location', 'southeast')