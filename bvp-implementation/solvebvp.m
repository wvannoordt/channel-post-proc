clear
clc
close all
Pr = 0.72;
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

Tcs   = csvread('data/bal-cs-T.csv');
Tcs   = Tcs(:,2);

Tref   = 100.0;
muref  = 3e-4;
mucs   = muref*(Tcs/Tref).^(0.76);

R = 287.15;
gamma = 1.4;
forc = 20.75/0.498742216e-02;

rhs = 0*ycs;
rhs = rhs + diff_prof(rhocs.*rhovkcs, ycs);
rhs = rhs + (gamma*R/(gamma-1))*diff_prof(rhocs.*vppTppcs, ycs);
rhs = rhs - diff_prof(utau0cs, ycs);

rhs = rhs - diff_prof(vtau1cs, ycs);
rhs = rhs - diff_prof(wtau2cs, ycs);
rhs = rhs - forc*rhocs.*utildcs;

n = length(ycs);
lhs = zeros(n,n);
lhs(1,1) = 1.0;
lhs(n,n) = 1.0;
rhs(1) = Tref;
rhs(n) = Tref;
for j=2:n-1
    mul = 0.5*(mucs(j)+mucs(j-1));
    mur = 0.5*(mucs(j)+mucs(j+1));
    y0 = ycs(j-1);
    y1 = ycs(j);
    y2 = ycs(j+1);
    a0 = (mul/Pr)*(1/(y1-y0))*(2/(y2-y0));
    a1 = (mur/Pr)*(1/(y2-y1))*(2/(y2-y0));
    lhs(j,j-1) = a0;
    lhs(j,j)   = -a0-a1;
    lhs(j,j+1) = a1;
end
TT = lhs\rhs;
figure
plot(ycs, TT/max(TT));
hold on
plot(ycs, Tcs/max(Tcs));
legend('Sol', 'DNS');