function r2 = recover(phi_mod, phi_mod2, f)


rdot  = 3*10^8/(4*pi)*(phi_mod -phi_mod2)/(0.1);%lambda/(4*pi)*phi_mod*1/T;

r2 = rdot - lsqnonlin(@(xx)getMulPath(rdot, xx, 0.2, f, phi_mod), 0.5);

end