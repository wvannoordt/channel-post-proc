clear
clc
close all

gamma = 1.4;
R = 287.15;
Twall = 100.0;

fs = 48;
lw = 3;

lstr = {};
lstr{end+1} = 'WMLES-M1.5';
lstr{end+1} = 'WMLES-M3.5';
lstr{end+1} = 'WMLES-M6.0';
lstr{end+1} = 'DNS-M1.5';
lstr{end+1} = 'DNS-M3.5';
lstr{end+1} = 'DNS-M6.0';

m0 = get_chan_sol('champs/m0/m1p5.csv');
m1 = get_chan_sol('champs/m1/m3p5.csv');
m2 = get_chan_sol('champs/m2/m6p0.csv');

[csT0, csY0] = csread('purdue/m0/cs-T.csv');
[csU0, csY0] = csread('purdue/m0/cs-u.csv');
[csruu0, csY0] = csread('purdue/m0/cs-upp.csv');
[csrvv0, csY0] = csread('purdue/m0/cs-vpp.csv');
[csrww0, csY0] = csread('purdue/m0/cs-wpp.csv');
[csuTpp0, csY0] = csread('purdue/m0/cs-uTpp.csv');
[csvTpp0, csY0] = csread('purdue/m0/cs-vTpp.csv');
[cswTpp0, csY0] = csread('purdue/m0/cs-wTpp.csv');

[csT1, csY1] = csread('purdue/m1/cs-T.csv');
[csU1, csY1] = csread('purdue/m1/cs-u.csv');
[csruu1, csY1] = csread('purdue/m1/cs-upp.csv');
[csrvv1, csY1] = csread('purdue/m1/cs-vpp.csv');
[csrww1, csY1] = csread('purdue/m1/cs-wpp.csv');
[csuTpp1, csY1] = csread('purdue/m1/cs-uTpp.csv');
[csvTpp1, csY1] = csread('purdue/m1/cs-vTpp.csv');
[cswTpp1, csY1] = csread('purdue/m1/cs-wTpp.csv');

[csT2, csY2] = csread('purdue/m2/cs-T.csv');
[csU2, csY2] = csread('purdue/m2/cs-u.csv');
[csruu2, csY2] = csread('purdue/m2/cs-upp.csv');
[csrvv2, csY2] = csread('purdue/m2/cs-vpp.csv');
[csrww2, csY2] = csread('purdue/m2/cs-wpp.csv');
[csuTpp2, csY2] = csread('purdue/m2/cs-uTpp.csv');
[csvTpp2, csY2] = csread('purdue/m2/cs-vTpp.csv');
[cswTpp2, csY2] = csread('purdue/m2/cs-wTpp.csv');

wsize = [100 100 1800 1500]
lloc  = 'southeastoutside';

figure('Position', wsize)
semilogx(m0.y, m0.T, 'linewidth', lw);
hold on
semilogx(m1.y, m1.T, 'linewidth', lw);
semilogx(m2.y, m2.T, 'linewidth', lw);
semilogx(csY0, csT0*Twall, 'linewidth', lw);
semilogx(csY1, csT1*Twall, 'linewidth', lw);
semilogx(csY2, csT2, 'linewidth', lw);
title('T', 'fontsize', fs)
set(gca, 'fontsize', fs)
h = legend(lstr);
set(h, 'fontsize', fs)
set(h, 'location', lloc)
saveas(gcf, 'T.png')

figure('Position', wsize)
semilogx(m0.y, m0.u, 'linewidth', lw);
hold on
semilogx(m1.y, m1.u, 'linewidth', lw);
semilogx(m2.y, m2.u, 'linewidth', lw);
semilogx(csY0, csU0*sqrt(gamma*R*Twall), 'linewidth', lw);
semilogx(csY1, csU1*sqrt(gamma*R*Twall), 'linewidth', lw);
semilogx(csY2, csU2, 'linewidth', lw);
title('U', 'fontsize', fs)
set(gca, 'fontsize', fs)
h = legend(lstr);
set(h, 'fontsize', fs)
set(h, 'location', lloc)
saveas(gcf, 'U.png')

figure('Position', wsize)
semilogx(m0.y, m0.ruu, 'linewidth', lw);
hold on
semilogx(m1.y, m1.ruu, 'linewidth', lw);
semilogx(m2.y, m2.ruu, 'linewidth', lw);
semilogx(csY0, csruu0*(sqrt(gamma*R*Twall))^2, 'linewidth', lw);
semilogx(csY1, csruu1*(sqrt(gamma*R*Twall))^2, 'linewidth', lw);
semilogx(csY2, csruu2, 'linewidth', lw);
title('Ruu', 'fontsize', fs)
set(gca, 'fontsize', fs)
h = legend(lstr);
set(h, 'fontsize', fs)
set(h, 'location', lloc)
saveas(gcf, 'ruu.png')

figure('Position', wsize)
semilogx(m0.y, m0.rvv, 'linewidth', lw);
hold on
semilogx(m1.y, m1.rvv, 'linewidth', lw);
semilogx(m2.y, m2.rvv, 'linewidth', lw);
semilogx(csY0, csrvv0*(sqrt(gamma*R*Twall))^2, 'linewidth', lw);
semilogx(csY1, csrvv1*(sqrt(gamma*R*Twall))^2, 'linewidth', lw);
semilogx(csY2, csrvv2, 'linewidth', lw);
title('Rvv', 'fontsize', fs)
set(gca, 'fontsize', fs)
h = legend(lstr);
set(h, 'fontsize', fs)
set(h, 'location', lloc)
saveas(gcf, 'rvv.png')

figure('Position', wsize)
semilogx(m0.y, m0.rww, 'linewidth', lw);
hold on
semilogx(m1.y, m1.rww, 'linewidth', lw);
semilogx(m2.y, m2.rww, 'linewidth', lw);
semilogx(csY0, csrww0*(sqrt(gamma*R*Twall))^2, 'linewidth', lw);
semilogx(csY1, csrww1*(sqrt(gamma*R*Twall))^2, 'linewidth', lw);
semilogx(csY2, csrww2, 'linewidth', lw);
title('Rww', 'fontsize', fs)
set(gca, 'fontsize', fs)
h = legend(lstr);
set(h, 'fontsize', fs)
set(h, 'location', lloc)
saveas(gcf, 'rww.png')

figure('Position', wsize)
semilogx(m0.y, m0.uTpp, 'linewidth', lw);
hold on
semilogx(m1.y, m1.uTpp, 'linewidth', lw);
semilogx(m2.y, m2.uTpp, 'linewidth', lw);
semilogx(csY0, csuTpp0*(sqrt(gamma*R*Twall))*Twall, 'linewidth', lw);
semilogx(csY1, csuTpp1*(sqrt(gamma*R*Twall))*Twall, 'linewidth', lw);
semilogx(csY2, csuTpp2, 'linewidth', lw);
title('uTpp', 'fontsize', fs)
set(gca, 'fontsize', fs)
h = legend(lstr);
set(h, 'fontsize', fs)
set(h, 'location', lloc)
saveas(gcf, 'uTpp.png')

figure('Position', wsize)
semilogx(m0.y, m0.vTpp, 'linewidth', lw);
hold on
semilogx(m1.y, m1.vTpp, 'linewidth', lw);
semilogx(m2.y, m2.vTpp, 'linewidth', lw);
semilogx(csY0, csvTpp0*(sqrt(gamma*R*Twall))*Twall, 'linewidth', lw);
semilogx(csY1, csvTpp1*(sqrt(gamma*R*Twall))*Twall, 'linewidth', lw);
semilogx(csY2, csvTpp2, 'linewidth', lw);
title('vTpp', 'fontsize', fs)
set(gca, 'fontsize', fs)
h = legend(lstr);
set(h, 'fontsize', fs)
set(h, 'location', lloc)
saveas(gcf, 'vTpp.png')

figure('Position', wsize)
semilogx(m0.y, m0.wTpp, 'linewidth', lw);
hold on
semilogx(m1.y, m1.wTpp, 'linewidth', lw);
semilogx(m2.y, m2.wTpp, 'linewidth', lw);
semilogx(csY0, cswTpp0*(sqrt(gamma*R*Twall))*Twall, 'linewidth', lw);
semilogx(csY1, cswTpp1*(sqrt(gamma*R*Twall))*Twall, 'linewidth', lw);
semilogx(csY2, cswTpp2, 'linewidth', lw);
title('wTpp', 'fontsize', fs)
set(gca, 'fontsize', fs)
h = legend(lstr);
set(h, 'fontsize', fs)
set(h, 'location', lloc)
saveas(gcf, 'wTpp.png')
