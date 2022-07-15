clear
clc
close all

ycs = csvread('data/bal-cs-qy.csv');
ycs = ycs(:,1);

Ncs = length(ycs);

qycs = csvread('data/bal-cs-qy.csv');
qycs = qycs(:,2);

utau0cs = csvread('data/bal-cs-utau10.csv');
utau0cs = utau0cs(:,2);

vtau1cs = csvread('data/bal-cs-vtau11.csv');
vtau1cs = vtau1cs(:,2);

wtau2cs = csvread('data/bal-cs-wtau12.csv');
wtau2cs = wtau2cs(:,2);

vppTppcs = csvread('data/bal-cs-vTpp.csv');
vppTppcs = vppTppcs(:,2);

rhovkcs = csvread('data/bal-cs-rhovk.csv');
rhovkcs = rhovkcs(:,2);

rhocs   = csvread('data/bal-cs-rho.csv');
rhocs   = rhocs(:,2);

utildcs = csvread('data/bal-cs-utild.csv');
utildcs = utildcs(:,2);

mucs   = 

Tbarcs = 

R = 287.15;
gamma = 1.4;
forc = 20.75/0.498742216e-02;

rhs = 0*ycs;
rhs = rhs + diff_prof(rhocs.*rhovkcs, ycs);
rhs = rhs + (gamma*R/(gamma-1))*diff_prof(rhocs.*vppTppcs, ycs);
rhs = rhs + diff_prof(utau0cs, ycs);
rhs = rhs + diff_prof(vtau1cs, ycs);
rhs = rhs + diff_prof(wtau2cs, ycs);
figure
plot(ycs, rhs);