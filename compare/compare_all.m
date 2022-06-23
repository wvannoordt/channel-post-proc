clear
clc
close all

oxdata = load('ox/data.csv');
oxdata = oxdata(1:floor(end/2),:);
y_ox = oxdata(:,1);
utild_ox = oxdata(:,9);
Ttild_ox = oxdata(:,12);
c10_ox   = oxdata(:,41);
a00_ox   = oxdata(:,32);
a01_ox   = oxdata(:,33);
a02_ox   = oxdata(:,34);
a10_ox   = oxdata(:,35);
a11_ox   = oxdata(:,36);
a12_ox   = oxdata(:,37);
a20_ox   = oxdata(:,38);
a21_ox   = oxdata(:,39);
a22_ox   = oxdata(:,40);

inner_Scale(oxdata, 9,  'purdue/cs-u.csv', 'U tilde', 'png/utild.png', 1);
centerline_norm(oxdata, 12, 'purdue/cs-T.csv', 'T tilde', 'png/Ttild.png', 2);
create_fig(oxdata, 41, 'purdue/cs-c10.csv', 'C10', 'png/c10.png', 3);
create_fig(oxdata, 32, 'purdue/cs-a00.csv', 'A00', 'png/a00.png', 4);
create_fig(oxdata, 33, 'purdue/cs-a01.csv', 'A01', 'png/a01.png', 4);
create_fig(oxdata, 34, 'purdue/cs-a02.csv', 'A02', 'png/a02.png', 4);
create_fig(oxdata, 35, 'purdue/cs-a10.csv', 'A10', 'png/a10.png', 4);
create_fig(oxdata, 36, 'purdue/cs-a11.csv', 'A11', 'png/a11.png', 4);
create_fig(oxdata, 37, 'purdue/cs-a12.csv', 'A12', 'png/a12.png', 4);
create_fig(oxdata, 35, 'purdue/cs-a20.csv', 'A20', 'png/a20.png', 4);
create_fig(oxdata, 36, 'purdue/cs-a21.csv', 'A21', 'png/a21.png', 4);
create_fig(oxdata, 37, 'purdue/cs-a22.csv', 'A22', 'png/a22.png', 4);
create_fig(oxdata, 55, 'purdue/cs-b00.csv', 'B00', 'png/b00.png', 5);
create_fig(oxdata, 56, 'purdue/cs-b01.csv', 'B01', 'png/b01.png', 5);
create_fig(oxdata, 57, 'purdue/cs-b02.csv', 'B02', 'png/b02.png', 5);
create_fig(oxdata, 17, 'purdue/cs-upp.csv', 'upp', 'png/upp.png', 6);
create_fig(oxdata, 18, 'purdue/cs-vpp.csv', 'vpp', 'png/vpp.png', 6);
create_fig(oxdata, 19, 'purdue/cs-wpp.csv', 'wpp', 'png/wpp.png', 6);
create_fig(oxdata, 29, 'purdue/cs-uTpp.csv', 'uTpp', 'png/uTpp.png', 7);
create_fig(oxdata, 30, 'purdue/cs-vTpp.csv', 'vTpp', 'png/vTpp.png', 7);
create_fig(oxdata, 31, 'purdue/cs-wTpp.csv', 'wTpp', 'png/wTpp.png', 7);