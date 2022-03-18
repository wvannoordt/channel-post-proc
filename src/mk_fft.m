clear
clc
close all

nt = 12;
nj = 9;
nk = 96;

ij = 9;

x = 1;
p = 2;
u = 3;
v = 4;
w = 5;
t = 6;
r = 7;

var = p;


tmp = dlmread('output_fftp/xpuvwtr_nt007_j008_k086.csv');
avgfft = 0*fft(tmp(:, var));

count = 0;

for it=1:nt
  %for ij=1:nj
    for ik=1:nk
      
      rit = it-1;
      rij = ij-1;
      rik = ik-1;
      
      filename = ['output_fftp/xpuvwtr_nt' num2str(rit, '%0.3d') '_j' num2str(rij, '%0.3d') '_k' num2str(rik, '%0.3d') '.csv'];

      data = dlmread(filename);
      d = data(:, var);
      dft = fft(d);
      avgfft = avgfft + dft;
      count = count + 1;
      
      
    end
  %end
end

avgfft = avgfft / count;
semilogx(log(abs(avgfft(2:(end/2)))), 'LineWidth', 2);
  