clc;
close all;

%% minimum
[m1, i1] = min(rclosses.CMRSEmV)
[m2, i2] = min(rclosses.CMRSEmW)

%% plot
figure(1);
plot(rclosses.CMRSEmV, '--ob', 'LineWidth', 1.5);
hold on;
plot(rclosses.CMRSEmW, '--or', 'LineWidth', 1.5);