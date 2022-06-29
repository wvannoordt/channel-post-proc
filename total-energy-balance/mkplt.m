clear
clc
close all

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
rhovkcs = rhovkcs(:,2);

rhovTcs = csvread('bal-cs-rhovT.csv');
rhovTcs = heatfac*rhovTcs(:,2);

rhocs   = csvread('bal-cs-rho.csv');
rhocs   = rhocs(:,2);

utildcs = csvread('bal-cs-utild.csv');
utildcs = utildcs(:,2);
figure
%plot(y, rhovk);
hold on
plot(y,     qy, 'linewidth', 3);
plot(ycs, qycs, 'linewidth', 3);
h = legend('WMLES', 'DNS');
title('q_y')

figure
hold on
plot(y, utau0, 'linewidth', 3)
plot(ycs, utau0cs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
title('u tauyx')

figure
hold on
plot(y, vtau1, 'linewidth', 3)
plot(ycs, vtau1cs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
title('v tauyy')

figure
hold on
plot(y, wtau2, 'linewidth', 3)
plot(ycs, wtau2cs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
title('w tauyz')

figure
hold on
plot(y, rhovT, 'linewidth', 3)
plot(ycs, rhovTcs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
title('\rho vT')

figure
hold on
plot(y, rho.*utild, 'linewidth', 3)
plot(ycs, rhocs.*utildcs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
title('\rho u')

figure
hold on
plot(y, rhovk, 'linewidth', 3)
plot(ycs, rhovkcs, 'linewidth', 3)
h = legend('WMLES', 'DNS');
title('\rho v k')