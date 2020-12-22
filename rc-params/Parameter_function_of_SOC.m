clear;
clc;
close all;

R0 = [3.74E-03 4.31E-04 4.31E-04 4.31E-04 4.29E-04 4.31E-04 4.30E-04 4.31E-04 4.31E-04 4.29E-04];
SOC = [100 90 80 70 60 50 40 30 20 10];
R1 = [5.58E-04 3.28E-04 1.47E-04 5.67E-04 4.32E-04 4.92E-04 7.55E-05 6.42E-04 9.11E-04 0.003296846];
R2 = [3.43E-04 4.11E-04 5.66E-04 1.37E-04 5.61E-04 6.63E-05 2.98E-04 8.15E-05 8.72E-04 0.002409846];
R3 = [2.50E-05 6.86E-05 3.31E-04 4.28E-04 1.68E-04 1.39E-04 4.65E-04 4.72E-04 2.30E-04 0.003940923];
TAU1 = [7.92E+00 34.50655625 1.04E+02 8.510638298 62.11180124 15.21606817 2.74E+02 12.01634223 62.89308176 62.81407035];
TAU2 = [4.61E+01 7.204610951 17.92757261 0.883392226 11.48369316 0.959692898 2.055498458 1.01010101 13.30494944 2.534854246];
TAU3 = [2.27E+01 1.14E+02 1.429796969 52.93806247 0.909918107 1.32E+02 29.00232019 63.93861893 1.469291801 16.10824742];

load rc3_cf_fmincon_params.mat

figure
plot(SOC,R0, '--o','Linewidth', 1.5);
hold on;
plot(SOC,r0_cf, '--or','Linewidth', 1.5);
xlabel('SOC [%]')
ylabel('R0[Ohms]')
title('R0 at different SOC')
legend('Relaxation', 'Fmincon')


figure
plot(SOC,R1, '--o','Linewidth', 1.5);
hold on;
plot(SOC,r1_cf, '--or','Linewidth', 1.5);
xlabel('SOC [%]')
ylabel('R1[Ohms]')
title('R1 at different SOC')
legend('Relaxation', 'Fmincon')


figure
plot(SOC,R2, '--o','Linewidth', 1.5);
hold on;
plot(SOC,r2_cf, '--or','Linewidth', 1.5);
xlabel('SOC [%]')
ylabel('R2[Ohms]')
title('R2 at different SOC')
legend('Relaxation', 'Fmincon')

figure
plot(SOC,R3, '--o','Linewidth', 1.5);
hold on;
plot(SOC,r3_cf, '--or','Linewidth', 1.5);
xlabel('SOC [%]')
ylabel('R3[Ohms]')
title('R3 at different SOC')
legend('Relaxation', 'Fmincon')


figure
plot(SOC,TAU1, '--o','Linewidth', 1.5);
hold on;
plot(SOC,r1_cf .* c1_cf, '--or','Linewidth', 1.5);
xlabel('SOC [%]')
ylabel('TAU1[sec]')
title('TAU1 at different SOC')
legend('Relaxation', 'Fmincon')

figure
plot(SOC,TAU2, '--o','Linewidth', 1.5);
hold on;
plot(SOC,r2_cf .* c2_cf, '--or','Linewidth', 1.5);
xlabel('SOC [%]')
ylabel('TAU2[sec]')
title('TAU2 at different SOC')
legend('Relaxation', 'Fmincon')

figure
plot(SOC,TAU3, '--o','Linewidth', 1.5);
hold on;
plot(SOC,r3_cf .* c3_cf, '--or','Linewidth', 1.5);
xlabel('SOC [%]')
ylabel('TAU3[sec]')
title('TAU3 at different SOC')
legend('Relaxation', 'Fmincon')
