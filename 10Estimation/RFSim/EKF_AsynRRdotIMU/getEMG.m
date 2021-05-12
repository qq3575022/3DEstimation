function E = getEMG(y,x)

HN = NaN(6+8,1);

HN(1:6+8)  = getHxkMG(x);

%HN
E = y - HN;

end