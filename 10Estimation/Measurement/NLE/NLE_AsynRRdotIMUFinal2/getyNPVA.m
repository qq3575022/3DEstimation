function [y, i, j, l, index1, len1] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, RFtime, time, i, j, l, m, factor)

% m = m
% i = i
% j = j
% l = l
% timeTt = time(m)
% accTt = acc_time(i)
% gyroTt = gyro_time(j)
% RF_time = RFtime(l)

%factor = 1000;

if i <= length(acc_time) && j <= length(gyro_time) && l <= length(RFtime) && abs(time(m) - acc_time(i)) == 0.0000 && abs(time(m) - gyro_time(j)) == 0.0000 && abs(time(m) - RFtime(l)) == 0.0000
  index = 1;
  len = 9+8;
  y = NaN(9+8,1);
  y(1:9+8) = [r_sim1(l), rdot_sim1(l), r_sim2(l), rdot_sim2(l), r_sim3(l), rdot_sim3(l), r_sim4(l), rdot_sim4(l), mag_data_x(j), gyro_data_x(j), mag_data_y(j), gyro_data_y(j), mag_data_z(j), gyro_data_z(j), factor*acc_data_x(i), factor*acc_data_y(i), factor*acc_data_z(i)];
  i = i + 1;
  j = j + 1;
  l = l + 1;
  
elseif j <= length(gyro_time) &&  l <= length(RFtime) && abs(time(m) - gyro_time(j)) == 0.0000 && abs(time(m) - RFtime(l)) == 0.0000
  index = 2;
  len = 6+8;
  y = NaN(6+8,1);
  y(1:6+8) = [r_sim1(l), rdot_sim1(l), r_sim2(l), rdot_sim2(l), r_sim3(l), rdot_sim3(l), r_sim4(l), rdot_sim4(l),mag_data_x(j), gyro_data_x(j), mag_data_y(j), gyro_data_y(j), mag_data_z(j), gyro_data_z(j)];
  j = j + 1;
  l = l + 1;
  
elseif i <= length(acc_time)  && l <= length(RFtime) && abs(time(m) - acc_time(i)) == 0.0000 &&  abs(time(m) - RFtime(l)) == 0.0000
  index = 3;
  len = 3+8;
  y = NaN(3+8,1);
  y(1:3+8) = [r_sim1(l), rdot_sim1(l), r_sim2(l), rdot_sim2(l), r_sim3(l), rdot_sim3(l), r_sim4(l), rdot_sim4(l),factor*acc_data_x(i), factor*acc_data_y(i), factor*acc_data_z(i)];
  i = i + 1;
  l = l + 1;
  
elseif i <= length(acc_time) &&  j <= length(gyro_time) && abs(time(m) - acc_time(i)) ==0.0000 &&  abs(time(m) - gyro_time(j)) == 0.0000
  index = 4;
  len = 9;
  y = NaN(9,1);
  y(1:9) = [mag_data_x(j), gyro_data_x(j), mag_data_y(j), gyro_data_y(j), mag_data_z(j), gyro_data_z(j), factor*acc_data_x(i), factor*acc_data_y(i), factor*acc_data_z(i)];
  i = i + 1;
  j = j + 1;
  
  
elseif  l <= length(RFtime) && abs(time(m) - RFtime(l)) == 0.0000
  index = 5;
  len = 8;
  y = NaN(8,1);
  y(1:8) = [r_sim1(l), rdot_sim1(l), r_sim2(l), rdot_sim2(l), r_sim3(l), rdot_sim3(l), r_sim4(l), rdot_sim4(l)];
  l = l + 1;
  
elseif i <= length(acc_time) && abs(time(m) - acc_time(i)) == 0.0000
  index = 6;
  len = 3;
  y = NaN(3,1);
  y(1:3) = [factor*acc_data_x(i), factor*acc_data_y(i), factor*acc_data_z(i)];
  i = i + 1;
  
elseif j <= length(gyro_time)  && abs(time(m) - gyro_time(j)) == 0.0000
  index = 7;
  len = 6;
  y = NaN(6,1);
  y(1:6) = [mag_data_x(j), gyro_data_x(j), mag_data_y(j), gyro_data_y(j), mag_data_z(j), gyro_data_z(j)];
  j = j + 1;
  
else
    y = nan;
    len = -1;
    index = -1;
  
end

len1 = len;
index1 = index;


end