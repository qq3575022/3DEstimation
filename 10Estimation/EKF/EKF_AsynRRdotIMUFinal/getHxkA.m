function H = getHxkA(x,factor)

H = NaN(3,1);


H(1,1) =  x(3)*cos(x(14))*cos(x(12))+x(6)*(cos(x(14))*sin(x(12))*sin(x(10)) - sin(x(14))*cos(x(10))) + x(9)*(cos(x(14))*sin(x(12))*cos(x(10))+sin(x(14))*sin(x(10)));
H(2,1) =  x(3)*sin(x(14))*cos(x(12))+x(6)*(sin(x(14))*sin(x(12))*sin(x(10)) + cos(x(14))*cos(x(10))) + x(9)*(sin(x(14))*sin(x(12))*cos(x(10))-cos(x(14))*sin(x(10)));
H(3,1) = -x(3)*sin(x(12))           +x(6)*cos(x(12))*sin(x(10)) + x(9)*cos(x(12))*cos(x(10));


H(1,1) = factor*H(1,1);
H(2,1) = factor*H(2,1);
H(3,1) = factor*H(3,1);
end