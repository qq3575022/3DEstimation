function H = getHxkMA(x)

H = NaN(6,1);

H(1,1) = x(10);
H(2,1) = x(12);
H(3,1) = x(14);

H(4,1) =  x(3)*cos(x(14))*cos(x(12))+x(6)*(cos(x(14))*sin(x(12))*sin(x(10)) - sin(x(14))*cos(x(10))) + x(9)*(cos(x(14))*sin(x(12))*cos(x(10))+sin(x(14))*sin(x(10)));
H(5,1) =  x(3)*sin(x(14))*cos(x(12))+x(6)*(sin(x(14))*sin(x(12))*sin(x(10)) + cos(x(14))*cos(x(10))) + x(9)*(sin(x(14))*sin(x(12))*cos(x(10))-cos(x(14))*sin(x(10)));
H(6,1) = -x(3)*sin(x(12))           +x(6)*cos(x(12))*sin(x(10)) + x(9)*cos(x(12))*cos(x(10));

    
end