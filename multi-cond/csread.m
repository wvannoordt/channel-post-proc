function [csout, yout] = csread(filename)
  data  = csvread(filename);
  n = size(data, 1)/2;
  yout  = data(1:n,1);
  csout = data(1:n,2);
end