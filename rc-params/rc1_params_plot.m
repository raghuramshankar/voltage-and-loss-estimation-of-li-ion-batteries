clear;
clc;
close all;

%% load data
load rc2_cf_fmincon_params.mat
r0 = r0_cf;
r1 = r1_cf;
c1 = c1_cf;
load rc2_cf_params.mat
soc = [100; 90; 80; 70; 60; 50; 40; 30; 20; 10];

%% plot
figure(1);
subplot(2, 2, 1)
plot(soc, r0, '-ob', 'LineWidth', 2);
hold on;
plot(soc, r0_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('R0')
title('R0')
legend('Fmincon', 'Relaxation')

subplot(2, 2, 2)
plot(soc, r1, '-ob', 'LineWidth', 2);
hold on;
plot(soc, r1_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('R1')
title('R1')
legend('Fmincon', 'Relaxation')

subplot(2, 2, 3)
plot(soc, c1, '-ob', 'LineWidth', 2);
hold on;
plot(soc, c1_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('C1')
title('C1')
legend('Fmincon', 'Relaxation')