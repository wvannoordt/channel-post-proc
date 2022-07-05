clear
clc
close all

ycs = csvread('data/bal-cs-qy.csv');
ycs = ycs(:,1);

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

R = 287.15;
gamma = 1.4;
forc = 20.75/0.498742216e-02;