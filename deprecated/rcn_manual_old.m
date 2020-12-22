clear;
clc;
close all;

%% load data
% load('pulse_2_c.mat');
load('ocv_final.mat');
load('pulse_1_c_1.mat');

%% voltage limitations
v_max = 4.4;                                    % cut-over voltage
v_min = 3.0;                                    % cut-off voltage

%% resample data
t1 = time(1);
t2 = time(end);
dt = 1;                                         % sample time
I = [1e-5; I];
t = (t1:dt:t2) - t1;

I = interp1(time, I, t1:dt:t2);
V_batt_time = interp1(time, V_batt_time, t1:dt:t2);
ocv = interp1(time, ocv, t1:dt:t2);
soc = interp1(time, soc, t1:dt:t2);
s = length(t);

%% variables
Q = 16.8;
eta = 1.0;                                  % coulombic efficiency

n = 3;                                      % number of rc links
i_r = zeros(n, 1);
v_c = zeros(n, 1);
v_t = zeros(n, 1);
flag = 1;

%% compute soc
ind_init = interp1(soc, 1:length(soc), 1.0, 'nearest');
v_t(1) = ocv(ind_init);
z0 = soc(ind_init);
z = z0 - dt/(Q*3600)*eta*cumsum(I(1:end));

% manual fit
r0 = 3.74e-3;
r1 = 0.56e-3;
c1 = 14.2e3;
r2 = 0.34e-3;
c2 = 0.135e6;
r3 = 0.025e-3;
c3 = 0.907e6;

r = [r1; r2; r3];
c = [c1; c2; c3];

%% compute terminal voltage
for j = 1:n
    f(j) = exp(-dt/(r(j)*c(j)));
end
for k = 1:s-1
    if flag == 1
        i_r(:, k+1) = diag(f)*i_r(:, k) + (1-f')*I(k);
        v_c(:, k) = i_r(:, k).*r;
        ind = interp1(soc, 1:length(soc), z(k), 'nearest');
        v_t(k+1) = ocv(ind) - sum(v_c(:, k)) - I(k).*r0;
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
v_t = reshape(v_t, size(V_batt_time));
rmse = 1000*sqrt(mean((V_batt_time - v_t).^2))     % [mV]
ocv(ind)

%% plot
figure(1);
hold on;
plot(t, v_t, 'r', 'LineWidth', 2)
hold on;
plot(t, V_batt_time, 'b', 'LineWidth', 2)
grid on;
xlabel('Time [s]')
ylabel('Terminal Voltage [V]')
title('RC')
legend('RC', 'Physical', 'location', 'southeast')