function H = getHxkA(x)


x1 = [0,    0,  0.756];%[2.6256, 0.0889,0.858];%[-0.2,-0.2, 0.4];
x2 = [2.29, 0,  1.07];%[ 0.9, 0.8, 0.2];
x3 = [2.29,2.52,  0.756];%[ 0.8,-0.3,-0.2];
x4 = [0, 2.52,  1.07];

H = NaN(11,1);

H(1,1) = sqrt((x(1)-x1(1))^2+(x(4)-x1(2))^2+(x(7)-x1(3))^2);
H(2,1) = ((x(1)-x1(1))*x(2) +(x(4)-x1(2))*x(5) +(x(7)-x1(3))*x(8))/sqrt((x(1)-x1(1))^2+(x(4)-x1(2))^2+(x(7)-x1(3))^2);

H(3,1) = sqrt((x(1)-x2(1))^2+(x(4)-x2(2))^2+(x(7)-x2(3))^2);
H(4,1) = ((x(1)-x2(1))*x(2) +(x(4)-x2(2))*x(5) +(x(7)-x2(3))*x(8))/sqrt((x(1)-x2(1))^2+(x(4)-x2(2))^2+(x(7)-x2(3))^2);

H(5,1) = sqrt((x(1)-x3(1))^2+(x(4)-x3(2))^2+(x(7)-x3(3))^2);
H(6,1) = ((x(1)-x3(1))*x(2) +(x(4)-x3(2))*x(5) +(x(7)-x3(3))*x(8))/sqrt((x(1)-x3(1))^2+(x(4)-x3(2))^2+(x(7)-x3(3))^2);

H(7,1) = sqrt((x(1)-x4(1))^2+(x(4)-x4(2))^2+(x(7)-x4(3))^2);
H(8,1) = ((x(1)-x4(1))*x(2) +(x(4)-x4(2))*x(5) +(x(7)-x4(3))*x(8))/sqrt((x(1)-x4(1))^2+(x(4)-x4(2))^2+(x(7)-x4(3))^2);

H(9,1) =  x(3)*cos(x(14))*cos(x(12))+x(6)*(cos(x(14))*sin(x(12))*sin(x(10)) - sin(x(14))*cos(x(10))) + x(9)*(cos(x(14))*sin(x(12))*cos(x(10))+sin(x(14))*sin(x(10)));
H(10,1) =  x(3)*sin(x(14))*cos(x(12))+x(6)*(sin(x(14))*sin(x(12))*sin(x(10)) + cos(x(14))*cos(x(10))) + x(9)*(sin(x(14))*sin(x(12))*cos(x(10))-cos(x(14))*sin(x(10)));
H(11,1) = -x(3)*sin(x(12))           +x(6)*cos(x(12))*sin(x(10)) + x(9)*cos(x(12))*cos(x(10));

   
end