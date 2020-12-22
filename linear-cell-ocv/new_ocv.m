clear;
clc;
close all;

%% graphite
ocv_g = linspace(2.78, 0.02, 3691);
soc_g = linspace(0.0, 0.98, 3691);
figure(1);
plot(soc_g, ocv_g, 'linewidth', 2);
xlabel('SOC [%]')
ylabel('OCV [V]')
title('Linear Graphite')

%% NMC
ocv_n = linspace(4.44, 2.688, 3691);
soc_n = linspace(0.0, 0.975, 3691);
figure(2);  
plot(soc_n, ocv_n, 'linewidth', 2);
xlabel('SOC [%]')
ylabel('OCV [V]')
title('Linear NMC')

%% full cell
ocv_f = linspace(4.16, 3.0, 3691);
soc_f = linspace(0, 1, 3691);
figure(3);
plot(soc_f, ocv_f, 'linewidth', 2);
xlabel('SOC [%]')
ylabel('OCV [V]')
title('Linear Full Cell')