clear
clc
close all

gamma = 1.4;
R = 287.15;
Twall = 100.0;

m0 = get_chan_sol('champs/m0/m1p5.csv');
m1 = get_chan_sol('champs/m1/m3p5.csv');
m2 = get_chan_sol('champs/m2/m6p0.csv');

[csT0, csY0] = csread('purdue/m0/cs-T.csv');
[csU0, csY0] = csread('purdue/m0/cs-u.csv');

[csT1, csY1] = csread('purdue/m1/cs-T.csv');
[csU1, csY1] = csread('purdue/m1/cs-u.csv');

[csT2, csY2] = csread('purdue/m2/cs-T.csv');
[csU2, csY2] = csread('purdue/m2/cs-u.csv');

figure
semilogx(m0.y, m0.T);
hold on
semilogx(m1.y, m1.T);
semilogx(m2.y, m2.T);
semilogx(csY1, csT1*Twall);
semilogx(csY2, csT2);

title('T')

figure
semilogx(m0.y, m0.u);
hold on
semilogx(m1.y, m1.u);
semilogx(m2.y, m2.u);
semilogx(csY1, csU1*3.5*sqrt(gamma*R*Twall));
semilogx(csY2, csU2);
title('U')