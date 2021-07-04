function E = getEF(y,x, factor)

HN = NaN(4,1);

HN(1:4)  = getHxkF(x, factor);

%HN
E = y - HN;

end