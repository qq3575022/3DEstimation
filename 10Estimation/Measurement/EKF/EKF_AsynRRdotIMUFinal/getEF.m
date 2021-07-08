function E = getEF(y,x, factor)

HN = NaN(8,1);

HN(1:8)  = getHxkF(x, factor);

%HN
E = y - HN;

end