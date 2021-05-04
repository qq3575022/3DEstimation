clc, clear, close all

load('res.mat')
load('x.mat')

y = [res ones(length(res),1)];

A = eye(3,4);
A(1,1) = x(4); A(1,2) = x(1)*x(5);  A(1,3) = x(2)*x(6); A(1,4) = x(7);
               A(2,2) = x(5);       A(2,3) = x(3)*x(6); A(2,4) = x(8);
                                    A(3,3) = x(6);      A(3,4) = x(9);
                                    
yy = A * y';
yR = zeros(length(y),1);

for i = 1 : 1 : length(yR)
    yR(i) = sqrt(yy(1,i)^2 + yy(2,i)^2 + yy(3,i)^2) - 9.7932;
end

% error mean
mean(yR)
% error standard deviation
sqrt(var(yR))
