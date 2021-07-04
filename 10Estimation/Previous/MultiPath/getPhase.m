function H = getPhase(z1, xx, gamma, f)

lambda  = 3*10^8/f;

H = 2*pi*(z1 - xx)/lambda + atan(sin(2*pi*xx/lambda)/(z1/(gamma*(z1 - xx))+cos(2*pi*xx/lambda)));
% H(9,1) = x(3); %x(3)*cos(x(14))*cos(x(12))+x(6)*(cos(x(8))*sin(x(12))*sin(x(10)) - sin(x(14))*cos(x(10))) + x(9)*(cos(x(14))*sin(x(12))*cos(x(10))+sin(x(14))*sin(x(10)));
% H(10,1) = x(6); % x(3)*sin(x(14))*cos(x(12))+x(6)*(sin(x(8))*sin(x(12))*sin(x(10)) + cos(x(14))*cos(x(10))) + x(9)*(sin(x(14))*sin(x(12))*cos(x(10))-cos(x(14))*sin(x(10)));
% H(11,1) = x(9); %-x(3)*sin(x(12))           +x(6)*cos(x(12))*sin(x(10)) + x(9)*cos(x(12))*cos(x(10));
%     
end