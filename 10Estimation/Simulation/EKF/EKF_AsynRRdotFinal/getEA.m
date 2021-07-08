function E = getEA(y,x, factor)

HN = NaN(3,1);

HN(1:3)  = getHxkA(x, factor);

%HN
E = y - HN;

end