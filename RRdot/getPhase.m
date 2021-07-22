function [y1, y2, y3, y4] = getPhase(x1, x2, x3, x4, xx)

y1 = sqrt((x1(1) - xx(1))^2 + (x1(2) - xx(2))^2 +(x1(3) - xx(3))^2);
y2 = sqrt((x2(1) - xx(1))^2 + (x2(2) - xx(2))^2 +(x2(3) - xx(3))^2);
y3 = sqrt((x3(1) - xx(1))^2 + (x3(2) - xx(2))^2 +(x3(3) - xx(3))^2);
y4 = sqrt((x4(1) - xx(1))^2 + (x4(2) - xx(2))^2 +(x4(3) - xx(3))^2);


end