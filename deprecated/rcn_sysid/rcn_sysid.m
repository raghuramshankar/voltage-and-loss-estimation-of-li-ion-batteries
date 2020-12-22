clear;
clc;
close all;

%% load data
load('dis_c_20.mat');
load('pulse_10a.mat');

%% variables
s = size(time, 2);
n = 3;
dt = 5;

%% pulse discharge
% I = 1e-5*ones(floor(4*s/10), 1);
% I = [I; 10*ones(floor(s/10), 1)];
% I = [I; 10*ones(floor(5*s/10), 1)];

%% scrape voltage interval of interest
V_batt_time_rc = V_batt_time - 4.1364;
[m, i] = min(V_batt_time_rc);
V_batt_time_rc = - V_batt_time_rc(i) + V_batt_time_rc(i:end);
% V_batt_time_rc = V_batt_time_rc(i:end);
I_rc = [0; 10*ones(s-i, 1)];
time_rc = time(i:end);

i_r = zeros(n, 1);

%% sisosubid
% SISOsubid(V_batt_time_rc, I_rc, 2)

%% nonlinear least squares
x_0 = [1e-3; 1e-3; 1e3; 1e-3; 1e3; 1e-3; 1e3];
% r0 = 3.74e-3;
% r1 = 0.56e-3;
% c1 = 367.98;
% r2 = 0.34e-3;
% c2 = 63.16e3;
% r3 = 0.025e-3;
% c3 = 0.173e6;
% x_0 = [r0, r1, c1, r2, c2, r3, c3];
x_lsq = lsqnonlin(@(x) lsqnonlin_opt(x, V_batt_time_rc, I_rc, dt, s, i_r, n, i), x_0)

%% simulate rc
r0 = x_lsq(1);
for k = 1:n
    r(k) = x_lsq(2*k);
    c(k) = x_lsq(2*k+1);
end
for j = 1:n
    f(j) = exp(-dt/(r(j).*c(j)));
end
for k = 1:s-i+1
    i_r(:, k+1) = diag(f)*i_r(:, k) + (ones(n, 1)-f')*I_rc(k);
    v_c(:, k+1) = i_r(:,k).*r';
    v_t(k) = sum(v_c(:, k)) + I_rc(k).*r0;
end

rmse = 1000*sqrt(mean((V_batt_time_rc - v_t').^2))     % [mV]

%% plot
figure(1);
plot(time_rc, v_t, 'r', 'LineWidth', 2)
hold on;
plot(time_rc, V_batt_time_rc, 'b', 'LineWidth', 2)
xlabel('Time [s]')
ylabel('Terminal Voltage [V]')
title('RC1')
legend('RC1', 'Physical', 'location', 'southeast') 