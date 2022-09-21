function [s] = get_chan_sol(filename)
  data = csvread(filename);
  n = size(data,1)/2;
  s.y    = data(1:n, 1)+1.0;
  s.u    = data(1:n, 9);
  s.T    = data(1:n, 12);
  s.ruu  = data(1:n, 17);
  s.rvv  = data(1:n, 18);
  s.rww  = data(1:n, 19);
  s.uTpp = data(1:n, 29);
  s.vTpp = data(1:n, 30);
  s.wTpp = data(1:n, 31);
end