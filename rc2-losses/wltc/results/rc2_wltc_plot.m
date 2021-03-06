clear;
clc;
% close all;

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

%% load data 0.5c
load rc2_1_c.mat
v_t_manual = v_t;
loss_rc_manual = loss_rc;
load rc2_1_c_fmincon.mat
v_t_fmincon = v_t;
loss_rc_fmincon = loss_rc;
load v_t_zero

%% volt
figure
plot(time, v_t_zero, '-g', 'LineWidth', 1);
hold on;
plot(time, v_t_manual, '-c', 'LineWidth', 1);
hold on;
plot(time, v_t_fmincon, '-r', 'LineWidth', 1);
hold on;
plot(time, V_batt_time, '-b', 'LineWidth', 1);
hold on;
xlabel('Time [s]')
ylabel('Voltage [V]')
legend('OCV', 'Relaxation 1C', 'Full Pulse 1C', 'Physics-based', 'location', 'southwest')

%% loss
figure
plot(time, loss_phy, '-b', 'LineWidth', 1);
hold on;
plot(time, loss_rc_manual, '-c', 'LineWidth', 1);
hold on;
plot(time, loss_rc_fmincon, '-r', 'LineWidth', 1);
hold on;
xlabel('Time [s]')
ylabel('Loss [W]')
legend('Physics-based', 'Relaxation 1C', 'Full Pulse 1C', 'location', 'northwest')