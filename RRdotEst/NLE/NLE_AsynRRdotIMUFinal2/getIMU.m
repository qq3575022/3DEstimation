function [angle1, gyrox, gyroy, gyroz, accx, accy, accz, accT, gyroT, time] = getIMU()

data=readtable('2.csv','Delimiter', ',');  g=9.7953;   
angle = load('angle.mat').angle;%load('tdgyro.mat');

%
% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

time   = table2array(data(:,1));               
data_x = table2array(data(:,3));  data_y = table2array(data(:,5));  data_z = table2array(data(:,4));

acc_time   = time(find_acc_data);  acc_time   = (acc_time - acc_time(1))/1000000000;                
gyro_time  = time(find_gyro_data); gyro_time  = (gyro_time - gyro_time(1))/1000000000;
mag_time   = time(find_mag_data);  mag_time   = (mag_time - mag_time(1))/1000000000;

%%
% start time
xS  = find(abs(acc_time-107.99)<0.002);  xS = xS(1);    xE = find(abs(acc_time-111.984)<0.001);   xE = xE(1);
xxS = find(abs(gyro_time-107.99)<0.002); xxS = xxS(1); xxE = find(abs(gyro_time-111.984)<0.001); xxE = xxE(1);
xSS = find(abs(mag_time-107.99)<0.012);  xSS = xSS(1); xEE = find(abs(mag_time-111.984)<0.009);  xEE = xEE(1);
%time = unique(sort([accT; gyroT; magT]),'rows');

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x = medfilt1(data_x(find_acc_data), 100);
acc_data_y = -medfilt1(data_y(find_acc_data), 100);
acc_data_z = medfilt1(data_z(find_acc_data), 100)-g;

% AXY = [1, 0.0623, 0.0055; 0, 1, -0.0041; 0, 0, 1]*[0.9944, 0, 0; 0, 0.9999, 0; 0, 0, 0.9880]*([acc_data_x';acc_data_z';-acc_data_y'] + [0.1739; 0.0071; -0.2999]);
% 
% acc_data_x = AXY(1,:); acc_data_y = -AXY(3,:); acc_data_z = AXY(2,:);

%
gyro_data_x = medfilt1(data_x(find_gyro_data),100); gyro_data_y = medfilt1(data_y(find_gyro_data),100); gyro_data_z = medfilt1(data_z(find_gyro_data),100);
mag_data_x  = data_x(find_mag_data);   mag_data_y = data_y(find_mag_data);   mag_data_z = data_z(find_mag_data);

accx = acc_data_x(xS+35:xE+35) - mean(acc_data_x(xS+35:xE+35));% - 0.002*(acc_time(xS:xE) - acc_time(xS));
accy = acc_data_y(xS+35:xE+35) - mean(acc_data_y(xS+35:xE+35));% - 0.011*(acc_time(xS:xE) - acc_time(xS));
accz = acc_data_z(xS+35:xE+35) - mean(acc_data_z(xS+35:xE+35));
accT = acc_time(xS:xE);

gyrox = gyro_data_x(xxS+35:xxE+35) - mean(gyro_data_x(xxS+35:xxE+35));gyroy = gyro_data_y(xxS+35:xxE+35) - mean(gyro_data_y(xxS+35:xxE+35));gyroz = gyro_data_z(xxS+35:xxE+35) - mean(gyro_data_z(xxS+35:xxE+35));
gyroT = gyro_time(xxS:xxE);

magx = mag_data_x(xSS:xEE) - mean(mag_data_x(xSS:xEE));magy = mag_data_y(xSS:xEE) - mean(mag_data_y(xSS:xEE));magz = mag_data_z(xSS:xEE) - mean(mag_data_z(xSS:xEE));

angle1 = NaN(3, xxE - xxS + 1);

angle1(1,:) = angle(1,xxS:xxE) - mean(angle(1,xxS:xxE)); 
angle1(2,:) = angle(2,xxS:xxE) - mean(angle(2,xxS:xxE)); 
angle1(3,:) = angle(3,xxS:xxE) - mean(angle(3,xxS:xxE)); 

magT = mag_time(xSS:xEE);

% ========== Start End Time Acceleration ========== 
%%
time = unique(sort([accT; gyroT]),'rows');
%%
[PP1, VV1, AA1] = groundtruth1Dx2(time - time(1)); PP1 = PP1 - PP1(1);
[PP2, VV2, AA2] = groundtruth1Dy2(time - time(1));
[PP3, VV3, AA3] = groundtruth1Dz2(time - time(1));
PP2 = PP2 - PP2(1); PP3 = PP3 - PP3(1);

%%
coord3 = zeros(3, length(time));
% x
coord3(1,:) = PP1+1.03;
%y;
coord3(2,:) = PP2+1.31;
coord3(3,:) = PP3+1.03;


end