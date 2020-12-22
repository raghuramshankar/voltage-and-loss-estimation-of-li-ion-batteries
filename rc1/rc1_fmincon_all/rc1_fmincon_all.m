clear;
clc;
close all;

%% add paths
folder = fileparts(mfilename('fullpath')); 
addpath(genpath(folder));

%% run and save
tic
% rc1_fmincon_0;
% save('rc1_fmincon_0');
rc1_fmincon_1;
save('rc1_fmincon_1');
rc1_fmincon_2;
save('rc1_fmincon_2');
rc1_fmincon_3;
save('rc1_fmincon_3');
rc1_fmincon_3;
save('rc1_fmincon_3');
rc1_fmincon_4;
save('rc1_fmincon_4');
rc1_fmincon_5;
save('rc1_fmincon_5');
rc1_fmincon_6;
save('rc1_fmincon_6');
rc1_fmincon_7;
save('rc1_fmincon_7');
rc1_fmincon_8;
save('rc1_fmincon_8');
rc1_fmincon_9;
save('rc1_fmincon_9');
toc