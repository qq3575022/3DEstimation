function E = getEA(y,x)

HN = NaN(8,1);

HN(1:8)  = getHxkA(x);

%HN
E = y - HN;

end