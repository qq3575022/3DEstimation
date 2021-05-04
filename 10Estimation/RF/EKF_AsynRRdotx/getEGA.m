function E = getEGA(y,x)

HN = NaN(6+8,1);

HN(1:6+8)  = getHxkGA(x);
%HN

E = y - HN;

end