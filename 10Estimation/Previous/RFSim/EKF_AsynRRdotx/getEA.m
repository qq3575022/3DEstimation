function E = getEA(y,x)

HN = NaN(3+8,1);

HN(1:3+8)  = getHxkA(x);

%HN
E = y - HN;

end