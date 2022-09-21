clear
clc
close all

gamma = 1.4;
R = 287.15;
fac = gamma*R/(gamma - 1.0);

data1 = csvread('cs-alpha_t.csv');
data2 = csvread('cs-mu_t.csv');

ycs   = data1(:,1);
atcs  = data1(:,2);
mutcs = data2(:,2);
prtcs = mutcs./atcs;

n = 100;
semilogx(ycs(1:n), prtcs(1:n));