function E = getEG(y,x, factor)

HN = NaN(6,1);

HN(1:6)  = getHxkG(x, factor);

%HN
E = y - HN;

end