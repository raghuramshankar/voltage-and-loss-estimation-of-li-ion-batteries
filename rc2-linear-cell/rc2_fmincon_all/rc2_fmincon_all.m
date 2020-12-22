clear;
clc;
close all;

%% add paths
folder = fileparts(mfilename('fullpath')); 
addpath(genpath(folder));

%% run and save
tic
rc2_fmincon_0;
save('rc2_fmincon_0');
rc2_fmincon_1;
save('rc2_fmincon_1');
rc2_fmincon_2;
save('rc2_fmincon_2');
rc2_fmincon_3;
save('rc2_fmincon_3');
rc2_fmincon_3;
save('rc2_fmincon_3');
rc2_fmincon_4;
save('rc2_fmincon_4');
rc2_fmincon_5;
save('rc2_fmincon_5');
rc2_fmincon_6;
save('rc2_fmincon_6');
rc2_fmincon_7;
save('rc2_fmincon_7');
rc2_fmincon_8;
save('rc2_fmincon_8');
rc2_fmincon_9;
save('rc2_fmincon_9');
toc