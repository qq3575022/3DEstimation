function E = getEMGA(y,x)

HN = NaN(9+8,1);

HN(1:9+8)  = getHxkMGA(x);

%HN
E = y - HN;

end