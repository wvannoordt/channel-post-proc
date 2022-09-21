clear
clc
close all

set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

axis_fs  = 28;
label_fs = 28;

plot_lw  = 3;
axis_lw  = 2;

R = 287.15;
gamma = 1.4;
heatfac = gamma*R/(gamma-1);
data = csvread('data.csv');
n = 0.5*size(data,1);
y     = data(1:n, 1) + 1;
rho   = data(1:n, 3);
utild = data(1:n, 9);
rhovk = data(1:n, 61);
rhovT = heatfac*data(1:n, 62);
qy    = data(1:n, 63);
utau0 = data(1:n, 58);
vtau1 = data(1:n, 59);
wtau2 = data(1:n, 60);
rhovk = data(1:n, 61);

ycs = csvread('bal-cs-qy.csv');
ycs = ycs(:,1);

qycs = csvread('bal-cs-qy.csv');
qycs = qycs(:,2);

utau0cs = csvread('bal-cs-utau10.csv');
utau0cs = utau0cs(:,2);

vtau1cs = csvread('bal-cs-vtau11.csv');
vtau1cs = vtau1cs(:,2);

wtau2cs = csvread('bal-cs-wtau12.csv');
wtau2cs = wtau2cs(:,2);

rhovTcs = csvread('bal-cs-rhovT.csv');
rhovTcs = rhovTcs(:,2);

rhovkcs = csvread('bal-cs-rhovk.csv');
rhovkcs = 0.5*rhovkcs(:,2);

rhovTcs = csvread('bal-cs-rhovT.csv');
rhovTcs = heatfac*rhovTcs(:,2);

rhocs   = csvread('bal-cs-rho.csv');
rhocs   = rhocs(:,2);

utildcs = csvread('bal-cs-utild.csv');
utildcs = utildcs(:,2);
figure('Position', [500 500 1500 600])
hold on
plot(y,     qy, 'linewidth', 3);
plot(ycs, qycs, 'linewidth', 3);
h = legend('WMLES', 'DNS');
set(h, 'location', 'southeast')
set(h, 'fontsize', label_fs)
xlabel('$y/\delta$')
set(gca, 'linewidth', axis_lw)
set(gca, 'fontsize', axis_fs)
title('$q_y$', 'fontsize', label_fs)
saveas(gcf, 'qy.png')

figure('Position', [500 500 1500 600])
hold on
plot(y, utau0, 'linewidth', 3)
plot(ycs, utau0cs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
set(h, 'location', 'southeast')
set(h, 'fontsize', label_fs)
xlabel('$y/\delta$')
set(gca, 'linewidth', axis_lw)
set(gca, 'fontsize', axis_fs)
title('$u \tau_{yx}$', 'fontsize', label_fs)
saveas(gcf, 'utauyx.png')

% figure('Position', [500 500 1500 600])
% hold on
% plot(y, vtau1, 'linewidth', 3)
% plot(ycs, vtau1cs, 'linewidth', 3)
% h = legend('WMLES', 'DNS');
% set(h, 'location', 'southeast')
% set(h, 'fontsize', label_fs)
% title('$v \tau_{yy}$', 'fontsize', label_fs)

% figure('Position', [500 500 1500 600])
% hold on
% plot(y, wtau2, 'linewidth', 3)
% plot(ycs, wtau2cs, 'linewidth', 3)
% h = legend('WMLES', 'DNS');
% set(h, 'location', 'southeast')
% set(h, 'fontsize', label_fs)
% title('$w \tau_{yz}$', 'fontsize', label_fs)

figure('Position', [500 500 1500 600])
hold on
plot(y, rhovT, 'linewidth', 3)
plot(ycs, rhovTcs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
set(h, 'location', 'southeast')
set(h, 'fontsize', label_fs)
xlabel('$y/\delta$')
set(gca, 'linewidth', axis_lw)
set(gca, 'fontsize', axis_fs)
title('$\gamma R \overline{\rho vT}/ (\gamma - 1)$', 'fontsize', label_fs)
saveas(gcf, 'rhovT.png')

figure('Position', [500 500 1500 600])
hold on
plot(y, rho.*utild, 'linewidth', 3)
plot(ycs, rhocs.*utildcs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
set(h, 'location', 'southeast')
set(h, 'fontsize', label_fs)
xlabel('$y/\delta$')
set(gca, 'linewidth', axis_lw)
set(gca, 'fontsize', axis_fs)
title('$\overline{\rho u}$', 'fontsize', label_fs)
saveas(gcf, 'rhou.png')

figure('Position', [500 500 1500 600])
hold on
plot(y, rhovk, 'linewidth', 3)
plot(ycs, rhovkcs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
set(h, 'location', 'southeast')
set(h, 'fontsize', label_fs)
xlabel('$y/\delta$')
set(gca, 'linewidth', axis_lw)
set(gca, 'fontsize', axis_fs)
title('$\overline{\rho v k}$', 'fontsize', label_fs)
saveas(gcf, 'rhovk.png')

figure('Position', [500 500 1500 600])
hold on
plot(y, rho, 'linewidth', 3)
plot(ycs, rhocs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
set(h, 'location', 'southeast')
set(h, 'fontsize', label_fs)
xlabel('$y/\delta$')
set(gca, 'linewidth', axis_lw)
set(gca, 'fontsize', axis_fs)
title('$\rho$', 'fontsize', label_fs)
saveas(gcf, 'rho.png')