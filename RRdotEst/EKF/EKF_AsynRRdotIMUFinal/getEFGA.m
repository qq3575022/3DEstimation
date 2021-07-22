function E = getEFGA(y,x, factor)

HN = NaN(9+8,1);

HN(1:9+8)  = getHxkFGA(x, factor);

%HN
E = y - HN;

end