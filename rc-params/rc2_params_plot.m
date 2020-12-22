clear;
clc;
close all;

%% load data
load rc2_cf_fmincon_params.mat
r0 = r0_cf * 1e3;
r1 = r1_cf * 1e3;
r2 = r2_cf * 1e3;
c1 = c1_cf * 1e-3;
c2 = c2_cf * 1e-3;
load rc2_cf_params.mat
r0_cf = r0_cf * 1e3;
r1_cf = r1_cf * 1e3;
r2_cf = r2_cf * 1e3;
c1_cf = c1_cf * 1e-3;
c2_cf = c2_cf * 1e-3;
soc = [100; 90; 80; 70; 60; 50; 40; 30; 20; 10];

%% plot
figure(1);
subplot(2, 3, 1)
plot(soc, r0, '-ob', 'LineWidth', 2);
hold on;
plot(soc, r0_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('R0 [mOhm]')
title('R0')
legend('Fmincon', 'Relaxation')

subplot(2, 3, 2)
plot(soc, r1, '-ob', 'LineWidth', 2);
hold on;
plot(soc, r1_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('R1 [mOhm]')
title('R1')
legend('Fmincon', 'Relaxation')

subplot(2, 3, 3)
plot(soc, r2, '-ob', 'LineWidth', 2);
hold on;
plot(soc, r2_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('R2 [mOhm]')
title('R2')
legend('Fmincon', 'Relaxation')

subplot(2, 3, 4)
plot(soc, c1, '-ob', 'LineWidth', 2);
hold on;
plot(soc, c1_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('C1 [kF]')
title('C1')
legend('Fmincon', 'Relaxation')

subplot(2, 3, 5)
plot(soc, c2, '-ob', 'LineWidth', 2);
hold on;
plot(soc, c2_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('C2 [kF]')
title('C2')
legend('Fmincon', 'Relaxation')