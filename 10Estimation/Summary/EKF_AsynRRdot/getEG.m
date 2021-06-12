function E = getEG(y,x)

HN = NaN(3+8,1);

HN(1:3+8)  = getHxkG(x);

%HN
E = y - HN;

end