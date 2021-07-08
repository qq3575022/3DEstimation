function E = getEFG(y,x, factor)

HN = NaN(6+8,1);

HN(1:6+8)  = getHxkFG(x, factor);

%HN
E = y - HN;

end