clear;
clc;
close all;

%% standard figure
hight = 3.25; width = 3.25;
top = 0.2; bottom = 0.5; left = 0.5; right = 0.2;
set(0,'defaultFigureUnits','inches');
set(0,'defaultFigurePosition',[5 3 width hight]);
set(0,'defaultAxesLineWidth',0.5);
set(0,'defaultAxesYGrid','on');
set(0,'defaultAxesXGrid','on');
set(0,'defaultAxesFontName','Times');
set(0,'defaultAxesFontSize',10);
set(0,'defaultTextFontName','Times');
set(0,'defaultTextFontSize',10);
set(0,'defaultAxesUnits','normalized');
set(0,'defaultAxesPosition',[left/width bottom/hight (width-left-right)/width  (hight-bottom-top)/hight])
set(0,'defaultLineLineWidth',1.5);
grey = 0.5;
titlestart = -0.2;

%% load data
load('ocv_final.mat');
time_ocv = time;
% load('rc2_cf_0_5_c.mat');
% load('rc2_cf_fmincon_0_5_c.mat');
% load('rc2_cf_params.mat');
% load('rc2_cf_fmincon_params.mat');
% load('rc2_cf_5_c.mat');
% load('rc2_cf_fmincon_5_c.mat');
load('rc2_cf_10_c.mat');
% load('rc2_cf_fmincon_10_c.mat');
load('wltc.mat');
I = I + 1e-4;
I = I*1.95;

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

n = 2;                                      % number of rc links
i_r = zeros(n, 1);
v_c = zeros(n, 1);

%% compute soc
ind_init = interp1(soc, 1:length(soc), 1.0, 'nearest');
v_t(1) = ocv(ind_init);
v_t_zero(1) = ocv(ind_init);
z0 = soc(ind_init);
z = z0 - dt/(Q*3600)*eta*cumsum(I(1:end));

%% compute terminal voltage without rc
for k = 1:s-1
    ind = interp1(soc, 1:length(soc), z(k), 'nearest');
    v_t_zero(k+1) = ocv(ind);
end
v_t_zero = v_t_zero';

%% compute terminal voltage with rc
for k = 1:s-1
    p = 10 * (1 - z(k)) + 1;
%     r0(k) = interp1(r0_cf, p, 'linear');
%     r1(k) = interp1(r1_cf, p, 'linear');
%     r2(k) = interp1(r2_cf, p, 'linear');
%     c1(k) = interp1(c1_cf, p, 'linear');
%     c2(k) = interp1(c2_cf, p, 'linear');
   
    r0(k) = interp1(r0_cf, p, 'nearest');
    r1(k) = interp1(r1_cf, p, 'nearest');
    r2(k) = interp1(r2_cf, p, 'nearest');
    c1(k) = interp1(c1_cf, p, 'nearest');
    c2(k) = interp1(c2_cf, p, 'nearest');
    
    if p > (length(r0_cf) - 1)
        r0(k) = r0_cf(end - 1);
        r1(k) = r1_cf(end - 1);
        r2(k) = r2_cf(end - 1);
        c1(k) = c1_cf(end - 1);
        c2(k) = c2_cf(end - 1);
    end
    
    r = [r1(k); r2(k)];
    c = [c1(k); c2(k)];
    for j = 1:n
        f(j) = exp(-dt/(r(j)*c(j)));
    end
    i_r(:, k+1) = diag(f)*i_r(:, k) + (1-f')*I(k);
    v_c(:, k) = i_r(:, k).*r;
    ind = interp1(soc, 1:length(soc), z(k), 'nearest');
    v_t(k+1) = ocv(ind) - sum(v_c(:, k)) - I(k).*r0(k);
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
hold on;
plot(time, v_t, 'r', 'LineWidth', 2)
hold on;
plot(time, v_t_zero, 'g', 'LineWidth', 2)
hold on;
plot(time, V_batt_time, 'b', 'LineWidth', 2)
grid on;
xlabel('Time [s]')
ylabel('Terminal Voltage [V]')
% title('RC2 - WLTC')
legend('RC2', 'OCV', 'Physical', 'location', 'northeast')

figure(2);
plot(time, loss_rc, 'r', 'LineWidth', 2)
hold on;
plot(time, loss_phy, 'b', 'LineWidth', 2)
grid on;
xlabel('Time [s]')
ylabel('Loss [W]')
% title('Losses in RC2 and Physical Model - WLTC')
legend('RC2', 'Physical', 'location', 'northwest')