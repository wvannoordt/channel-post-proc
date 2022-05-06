clear
clc
close all

oxdata = load('ox/data.csv');
oxdata = oxdata(1:floor(end/2),:);
y_ox = oxdata(:,1);
utild_ox = oxdata(:,9);
Ttild_ox = oxdata(:,12);

utild_cs = load('purdue/cs-u.dat');
y_cs = utild_cs(:,1)-1;
utild_cs = utild_cs(:,2);

Ttild_cs = load('purdue/cs-T.dat');
Ttild_cs = Ttild_cs(:,2);

figure
hold on
plot(y_ox, utild_ox)
plot(y_cs, utild_cs)

figure
hold on
plot(y_ox, Ttild_ox)
plot(y_cs, Ttild_cs)
