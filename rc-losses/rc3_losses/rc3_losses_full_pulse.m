clear;
clc;
close all;

%% load data
load('ocv_final.mat');
time_ocv = time;
load('rc3_cf_params.mat');
load('full_pulse.mat');

%% voltage limitations
v_max = 4.4;                                    % cut-over voltage
v_min = 3.0;                                    % cut-off voltage

%% resample ocv and soc
t1 = time_ocv(1);
t2 = time_ocv(end);
dt = 1;                                         % sample time
t = (t1:dt:t2) - t1;

ocv = interp1(time_ocv, ocv, t1:dt:t2);
soc = interp1(time_ocv, soc, t1:dt:t2);

s = length(time);

%% variables
Q = 16.8;
eta = 1.0;                                  % coulombic efficiency

n = 3;                                      % number of rc links
i_r = zeros(n, 1);
v_c = zeros(n, 1);

%% compute soc
ind_init = interp1(soc, 1:length(soc), 1.0, 'nearest');
v_t(1) = ocv(ind_init);
v_t_zero(1) = ocv(ind_init);
z0 = soc(ind_init);
z = z0 - dt/(Q*3600)*eta*cumsum(I(1:end));

% % manual fit
% p = 1;
% r0 = r0_cf(p);
% r1 = r1_cf(p);
% r2 = r2_cf(p);
% r3 = r3_cf(p);
% c1 = c1_cf(p);
% c2 = c2_cf(p);
% c3 = c3_cf(p);

% r = [r1; r2; r3];
% c = [c1; c2; c3];

%% compute terminal voltage without rc
for k = 1:s-1
    ind = interp1(soc, 1:length(soc), z(k), 'nearest');
    v_t_zero(k+1) = ocv(ind);
end
v_t_zero = v_t_zero';

%% compute terminal voltage with rc
for k = 1:s-1
    p = 10 * (1 - z(k)) + 1;
%     r0 = interp1(r0_cf, p, 'linear');
%     r1 = interp1(r1_cf, p, 'linear');
%     r2 = interp1(r2_cf, p, 'linear');
%     r3 = interp1(r3_cf, p, 'linear');
%     c1 = interp1(c1_cf, p, 'linear');
%     c2 = interp1(c2_cf, p, 'linear');
%     c3 = interp1(c3_cf, p, 'linear');
    
    r0 = interp1(r0_cf, p, 'nearest');
    r1 = interp1(r1_cf, p, 'nearest');
    r2 = interp1(r2_cf, p, 'nearest');
    r3 = interp1(r3_cf, p, 'nearest');
    c1 = interp1(c1_cf, p, 'nearest');
    c2 = interp1(c2_cf, p, 'nearest');
    c3 = interp1(c3_cf, p, 'nearest');
    
    if p > 10
        r0 = r0_cf(10);
        r1 = r1_cf(10);
        r2 = r2_cf(10);
        r3 = r3_cf(10);
        c1 = c1_cf(10);
        c2 = c2_cf(10);
        c3 = c3_cf(10);
    end
    
    r = [r1; r2; r3];
    c = [c1; c2; c3];
    for j = 1:n
        f(j) = exp(-dt/(r(j)*c(j)));
    end
    i_r(:, k+1) = diag(f)*i_r(:, k) + (1-f')*I(k);
    v_c(:, k) = i_r(:, k).*r;
    ind = interp1(soc, 1:length(soc), z(k), 'nearest');
    v_t(k+1) = ocv(ind) - sum(v_c(:, k)) - I(k).*r0;
end

%% compute rms error
v_t = reshape(v_t, size(V_batt_time));
rmse_volt = 1000*sqrt(mean((V_batt_time - v_t).^2)) % [mV]

%% compute losses
op_rc = v_t_zero - v_t;
op_phy = v_t_zero - V_batt_time;

loss_rc = op_rc .* I;                               % [W]
loss_phy = op_phy .* I;                             % [W]

rmse_loss = 1000*sqrt(mean((loss_phy - loss_rc).^2))% [mW]

energy_rc = trapz(time, loss_rc)/3600               % [Wh]
energy_phy = trapz(time, loss_phy)/3600             % [Wh]

%% plot
figure(1);
plot(time, v_t, 'r', 'LineWidth', 2)
hold on;
plot(time, v_t_zero, 'g', 'LineWidth', 2)
hold on;
plot(time, V_batt_time, 'b', 'LineWidth', 2)
grid on;
xlabel('Time [s]')
ylabel('Terminal Voltage [V]')
title('RC - Full Pulse Discharge')
legend('RC', 'No RC', 'Physical', 'location', 'northeast')

figure(2);
plot(time, loss_rc, 'r', 'LineWidth', 2)
hold on;
plot(time, loss_phy, 'g', 'LineWidth', 2)
hold on;
xlabel('Time [s]')
ylabel('Loss [W]')
title('Losses in RC and Physical Model - Full Pulse')
legend('RC', 'Physical', 'location', 'northeast')