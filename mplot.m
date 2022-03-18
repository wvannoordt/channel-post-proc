clear
clc
close all

data = dlmread('output/data.csv');


y       = 1;
mu_bar  = 2;
rho_bar = 3;
u_bar   = 4;
v_bar   = 5;
w_bar   = 6;
T_bar   = 7;
P_bar   = 8;
u_tilde = 9;
v_tilde = 10;
w_tilde = 11;
T_tilde = 12;
up      = 13;
vp      = 14;
wp      = 15;
Tp      = 16;
upp     = 17;
vpp     = 18;
wpp     = 19;
Tpp     = 20;
  
figure
hold on
plot(data(:,y), data(:,Tpp));
plot(data(:,y), data(:,Tp));