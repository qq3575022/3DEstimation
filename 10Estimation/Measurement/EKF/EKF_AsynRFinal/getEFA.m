function E = getEFA(y,x, factor)

HN = NaN(11,1);

HN(1:11)  = getHxkFA(x, factor);

%HN
E = y - HN;

end