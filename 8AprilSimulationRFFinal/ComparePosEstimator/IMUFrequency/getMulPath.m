function E = getMulPath(r_sim1, r_sim2, r_sim3, r_sim4, x1, x2, x3, x4, xx)

%phi1 = mod(phi, 2*pi);

[y1, y2, y3, y4] = getPhase(x1, x2, x3, x4, xx);

%y1= mod(y, 2*pi);

E = (r_sim1 - y1).^2 + (r_sim2 - y2).^2 + (r_sim3 - y3).^2 + (r_sim4 - y4).^2;

end
