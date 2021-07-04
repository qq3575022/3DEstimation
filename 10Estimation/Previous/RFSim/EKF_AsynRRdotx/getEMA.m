function E = getEMA(y,x)

HN = NaN(6+8,1);

HN(1:6+8)  = getHxkMA(x);

%HN
E = y - HN;

end