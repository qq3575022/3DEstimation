
z1 = 5;
gamma = 1;
f = 5.8*10^9;
phi = 545.6808335;

x = lsqnonlin(@(xx)getMulPath(z1, xx, gamma, f, phi), 0.1);
