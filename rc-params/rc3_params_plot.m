clear;
clc;
close all;

%% load data
load rc3_cf_fmincon_params.mat
r1 = r1_cf;
r2 = r2_cf;

r3 = r3_cf;
c1 = c1_cf;
c2 = c2_cf;
c3 = c3_cf;
load rc3_cf_params.mat
soc = [100; 90; 80; 70; 60; 50; 40; 30; 20; 10];

%% plot
figure(1);
subplot(2, 3, 1)
plot(soc, r1, '-ob', 'LineWidth', 2);
hold on;
plot(soc, r1_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('R1')
title('R1')
legend('Fmincon', 'Relaxation')

subplot(2, 3, 2)
plot(soc, r2, '-ob', 'LineWidth', 2);
hold on;
plot(soc, r2_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('R2')
title('R2')
legend('Fmincon', 'Relaxation')

subplot(2, 3, 3)
plot(soc, r3, '-ob', 'LineWidth', 2);
hold on;
plot(soc, r3_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('R3')
title('R3')
legend('Fmincon', 'Relaxation')

subplot(2, 3, 4)
plot(soc, c1, '-ob', 'LineWidth', 2);
hold on;
plot(soc, c1_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('C1')
title('C1')
legend('Fmincon', 'Relaxation')

subplot(2, 3, 5)
plot(soc, c2, '-ob', 'LineWidth', 2);
hold on;
plot(soc, c2_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('C2')
title('C2')
legend('Fmincon', 'Relaxation')

subplot(2, 3, 6)
plot(soc, c3, '-ob', 'LineWidth', 2);
hold on;
plot(soc, c3_cf, '--or', 'LineWidth', 2);
xlabel('SOC [%]')
ylabel('C3')
title('C3')
legend('Fmincon', 'Relaxation')