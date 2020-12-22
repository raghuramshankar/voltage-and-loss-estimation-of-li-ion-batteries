clear;
clc;
close all;

%% load data
load('dis_c_20.mat');
ocv = V_batt_time;
load('pulse_10a.mat')

%% time and voltage limitations
v_max = 4.2;                                % cut-over voltage
v_min = 3;                                  % cut-off voltage
s = size(time, 2);

%% variables
dt = 5;                                     % sample time [s]
Q = 16.4;
eta = 1;                                    % coulombic efficiency

i_r = zeros(s-1, 1);
v_c = zeros(s-1, 1);
ind = zeros(s-1, 1);
v_t = zeros(s-1, 1);
v_t(1) = ocv(1);
soc = linspace(0.9, 0.1, s-1);
flag = 1;

%% pulse discharge
I = 1e-5*ones(floor(4*s/10), 1);
I = [I; 10*ones(floor(s/10), 1)];
I = [I; 1e-5*ones(floor(5*s/10), 1)];

%% compute soc
z0 = soc(1);
for k = 1:s-1
    z(k) = z0 - dt/(Q*3600)*eta*sum(I(1:k));
end

%% scrape voltage interval of interest
% V_batt_time_rc = V_batt_time - ocv(interp1(soc, 1:length(soc), z(end), 'nearest'));
% [m, i] = min(V_batt_time_rc);
% % V_batt_time_rc = - V_batt_time_rc(i) + V_batt_time_rc(i:end);
% V_batt_time_rc = V_batt_time_rc(i:end);
% time_rc = time(i:end);

%% fmincon
% constraints
x_0 = [0, 0, 0];                        % initial values [r0, r1, c1]
A = [];                                 % linear inequality constraints
B = [];
Aeq = [];                               % linear equality constraints
Beq = [];
lb = [0, 0, 0];                         % upper and lower bounds
ub = [1e8, 1e8, 1e8];

% run fmincon
fun = @(x) opt(x, V_batt_time, I, z, soc, dt, s, i_r, v_c, ind, v_t, ocv);
rc1_param = fmincon(fun,x_0,A,B,Aeq,Beq,lb,ub);
r0 = rc1_param(1)                           % [ohm]
r = rc1_param(2)                            % [ohm]
c = rc1_param(3)                            % [farad]

%% compute terminal voltage
f = exp(-dt/(r*c));
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
rmse = 1000*sqrt(mean((V_batt_time - v_t).^2))     % [mV]

%% plot
figure(1);
plot(time, v_t, 'r', 'LineWidth', 2)
hold on;
plot(time, V_batt_time, 'b', 'LineWidth', 2)
xlabel('Time [s]')
ylabel('Terminal Voltage [V]')
title('RC1 - fmincon fit')
legend('RC1', 'Physical', 'location', 'southeast')