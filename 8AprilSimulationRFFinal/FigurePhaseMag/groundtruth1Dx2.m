function [PP, VV, AA] = groundtruth1Dx2(td)  
%x 
p1 = 0;
p2 = 0.384;

ac = 0.128;
vc = 0.384/3;%0.128*sqrt(3);%0.384;

ta = vc/ac;
tv = (p2-p1)/vc - vc/ac;

% t1 = 0.157;
t2 = ta+tv+ta;

PP = NaN(1,length(td));
VV = NaN(1,length(td));
AA = NaN(1,length(td));

for n = 1:length(td)
    t = td(n);
    if t < ta
        p = p1+0.5*ac*t^2;
        v = ac*t;
        a = ac;
    elseif t < t2-ta
        p = 0.5*(p1+p2)+vc*(t-0.5*t2);
        v = vc;
        a = 0;
    elseif t < t2
        p = p2-0.5*ac*(t2-t)^2;
        v = -ac*(t-t2);
        a = -ac;
    else
        p = p2;
        v = 0;
        a = 0;
    end
    PP(n) = p;
    VV(n) = v;
    AA(n) = a;
end
end
    