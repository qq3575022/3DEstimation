function [y, index, len] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,  m)
% m = m
% i = i
% j = j
% k = k
% timeTt = time(m)
% accTt = acc_time(i)
% gyroTt = gyro_time(j)
% magTt = mag_time(k)

index = 1;
len = 8;
y = NaN(8,1);
y(1:8) = [r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4];

  
% if i <= length(acc_time) && j <= length(gyro_time) && k <= length(mag_time) && abs(time(m) - acc_time(i)) == 0 && abs(time(m) - gyro_time(j)) ==0 && abs(time(m) - mag_time(k)) ==0
%   index = 1;
%   len = 9+8;
%   y = NaN(9+8,1);
%   y(1:9+8) = [r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, mag_data_x(k), gyro_data_x(j), mag_data_y(k), gyro_data_y(j), mag_data_z(k), gyro_data_z(j), acc_data_x(i), acc_data_y(i), acc_data_z(i)];
%   i = i + 1;
%   j = j + 1;
%   k = k + 1;
%   
% elseif i <= length(acc_time) && j <= length(gyro_time) && abs(time(m) - acc_time(i)) == 0 && abs(time(m) - gyro_time(j)) == 0
%   index = 2;
%   len = 6+8;
%   y = NaN(6+8,1);
%   y(1:6+8) = [r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, gyro_data_x(j), gyro_data_y(j), gyro_data_z(j), acc_data_x(i), acc_data_y(i), acc_data_z(i)];
%   i = i + 1;
%   j = j + 1;
%   
% elseif j <= length(gyro_time) && k <= length(mag_time) && abs(time(m) - gyro_time(j)) == 0 && abs(time(m) - mag_time(k)) == 0
%   index = 3;
%   len = 6+8;
%   y = NaN(6+8,1);
%   y(1:6+8) = [r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, mag_data_x(k), gyro_data_x(j), mag_data_y(k), gyro_data_y(j), mag_data_z(k), gyro_data_z(j)];
%   j = j + 1;
%   k = k + 1;
%   
% elseif i <= length(acc_time) && k <= length(mag_time) && abs(time(m) - acc_time(i)) == 0 && abs(time(m) - mag_time(k)) == 0
%   index = 4;
%   len = 6+8;
%   y = NaN(6+8,1);
%   y(1:6+8) = [r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, mag_data_x(k), mag_data_y(k), mag_data_z(k), acc_data_x(i), acc_data_y(i), acc_data_z(i)];
%   i = i + 1;
%   k = k + 1;
%   
% elseif i <= length(acc_time) && abs(time(m) - acc_time(i)) == 0
%   index = 5;
%   len = 3+8;
%   y = NaN(3+8,1);
%   y(1:3+8) = [r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, acc_data_x(i), acc_data_y(i), acc_data_z(i)];
%   i = i + 1;
%   
% elseif j <= length(gyro_time) && abs(time(m) - gyro_time(j)) < 0.0001
%   index = 6;
%   len = 3+8;
%   y = NaN(3+8,1);
%   y(1:3+8) = [r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, gyro_data_x(j), gyro_data_y(j), gyro_data_z(j)];
%   j = j + 1;
%   
% elseif k <= length(mag_time) && abs(time(m) - mag_time(k)) < 0.0001
%   index = 7;
%   len = 3+8;
%   y = NaN(3+8,1);
%   y(1:3+8) = [r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, mag_data_x(k), mag_data_y(k), mag_data_z(k)];
%   k = k + 1;
%   
% else
%     y = nan;
%     len = -1;
%     index = -1;
%   
% end

end