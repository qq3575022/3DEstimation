function E = getEGA(y,x, factor)

HN = NaN(9,1);

HN(1:9)  = getHxkGA(x, factor);
%HN

E = y - HN;

end