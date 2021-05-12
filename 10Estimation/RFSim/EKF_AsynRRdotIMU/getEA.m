function E = getEA(y,x)

HN = NaN(11,1);

HN(1:11)  = getHxkA(x);

%HN
E = y - HN;

end