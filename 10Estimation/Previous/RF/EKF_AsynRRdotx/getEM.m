function E = getEM(y,x)

HN = NaN(3+8,1);

HN(1:3+8)  = getHxkM(x);

%HN
E = y - HN;

end