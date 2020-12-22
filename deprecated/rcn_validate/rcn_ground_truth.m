clear;
clc;
close all;

%% variables
s = 100;
n = 3;
dt = 1;

%% pulse discharge
I_rc = [0; 10*ones(s-1, 1)];
time_rc = linspace(0, s, s);

i_r = zeros(n, 1);

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

%% simulate rc
r0 = x_0(1);
for k = 1:n
    r(k) = x_0(2*k);
    c(k) = x_0(2*k+1);
end
for j = 1:n
    f(j) = exp(-dt/(r(j).*c(j)));
end
for k = 1:s-1
    i_r(:, k+1) = diag(f)*i_r(:, k) + (ones(n, 1)-f')*I_rc(k);
    v_c(:, k+1) = i_r(:,k).*r';
    v_t(k+1) = sum(v_c(:, k)) + I_rc(k).*r0;
end

%% plot
figure(1);
plot(time_rc, v_t, 'r', 'LineWidth', 2)
xlabel('Time [s]')
ylabel('Terminal Voltage [V]')
title('RC')
legend('RC', 'location', 'southeast') 