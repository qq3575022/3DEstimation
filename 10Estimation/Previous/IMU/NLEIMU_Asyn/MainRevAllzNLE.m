clc, clear, close all
data=readtable('2.csv','Delimiter', ',');  g=9.7953;
load('tdgyro.mat'); load('angle.mat');

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

time   = table2array(data(:,1));               
data_x = table2array(data(:,3));  data_y = table2array(data(:,5));  data_z = table2array(data(:,4));

acc_time   = time(find_acc_data);  acc_time   = (acc_time - acc_time(1))/1000000000;                
gyro_time  = time(find_gyro_data); gyro_time  = (gyro_time - gyro_time(1))/1000000000;
mag_time   = time(find_mag_data);  mag_time   = (mag_time - mag_time(1))/1000000000;

% start time
% 

zS = find(abs(acc_time-67.9776)<0.002); zS = zS(1);
zE = find(abs(acc_time-70.0537)<0.002); zE = zE(1);

zzS = find(abs(gyro_time-67.9776)<0.002); zzS = zzS(1);
zzE = find(abs(gyro_time-70.0537)<0.002); zzE = zzE(1);

zSS = find(abs(mag_time-67.9776)<0.009); zSS = zSS(1);
zEE = find(abs(mag_time-70.0537)<0.009); zEE = zEE(1);
% 

%time = unique(sort([accT; gyroT; magT]),'rows');

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x = data_x(find_acc_data);
acc_data_y = -data_y(find_acc_data);
acc_data_z = data_z(find_acc_data)-g;

AXY = [1, 0.0623, 0.0055; 0, 1, -0.0041; 0, 0, 1]*[0.9944, 0, 0; 0, 0.9999, 0; 0, 0, 0.9880]*([acc_data_x';acc_data_z';-acc_data_y'] + [0.1739; 0.0071; -0.2999]);

acc_data_x = AXY(1,:); acc_data_y = -AXY(3,:); acc_data_z = AXY(2,:);

gyro_data_x = data_x(find_gyro_data);gyro_data_y = data_y(find_gyro_data);gyro_data_z = data_z(find_gyro_data);
mag_data_x = data_x(find_mag_data);mag_data_y = data_y(find_mag_data);mag_data_z = data_z(find_mag_data);

% acc_data_x = medfilt1(acc_data_x,30);
% acc_data_y = medfilt1(acc_data_y,30);
% acc_data_z = medfilt1(acc_data_z,30);

accx = acc_data_x(zS:zE) - mean(acc_data_x(zS:zE));accy = acc_data_y(zS:zE) - mean(acc_data_y(zS:zE))+0.003;accz = acc_data_z(zS:zE) - mean(acc_data_z(zS:zE));
accT = acc_time(zS:zE);

gyrox = gyro_data_x(zzS:zzE) - mean(gyro_data_x(zzS:zzE));gyroy = gyro_data_y(zzS:zzE) - mean(gyro_data_y(zzS:zzE));gyroz = gyro_data_z(zzS:zzE) - mean(gyro_data_z(zzS:zzE));

angle(1,:) = angle(1,:) - mean(angle(1,:)); 
angle(2,:) = angle(2,:) - mean(angle(2,:)); 
angle(3,:) = angle(3,:) - mean(angle(3,:)); 


gyroT = gyro_time(zzS:zzE);

magx = mag_data_x(zSS:zEE) - mean(mag_data_x(zSS:zEE));magy = mag_data_y(zSS:zEE) - mean(mag_data_y(zSS:zEE));magz = mag_data_z(zSS:zEE) - mean(mag_data_z(zSS:zEE));
magT = mag_time(zSS:zEE);
% ========== Start End Time Acceleration ========== 

time = unique(sort([accT; gyroT; magT]),'rows');

[PP1, VV1, AA1] = groundtruth1Dz2(time - time(1));

% Start from 0
PP1 = PP1 - PP1(1);

% Measurement Vector [acc_x, acc_y, acc_z, orientation, angular velocity]

% State vector [xdotdot, ydotdot, orientation, angular velocity]
%x = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]'*ones(1,length(accT));
%o1 a1 o2 a2 o3 a3
N = 5; i = 1; j = 1; k = 1;
% Initiate x, e and temproary variable
x = [0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]*ones(1,length(time)-N);e = NaN(1,length(time)-N);

for m = 1:1:length(time)-N
    
    [y, i, j, k, index, len] = gety(angle(1,:), gyrox, angle(2,:), gyroy, angle(3,:), gyroz, accx, accy, accz, accT, gyroT, magT, time, i, j, k, m);

    x(:,m) = lsqnonlin(@(xx)getNLE(y, xx, N, index, len, time, m),[0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]);

    F = getF(time, m, m+1);
    x(:,m) = F*x(:,m);

    
end
%%
[pos, vel] = doubleInte(diff(time(1:end-N)), x(9,2:end));


figure
subplot(311), plot(time(2:end-N), pos, 'LineWidth', 3); hold on; plot(time, PP1, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','SouthEast'); title('Position Along z Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([67.9376, 70.0537]);
subplot(312), plot(time(2:end-N), vel, 'LineWidth', 2); hold on; plot(time, VV1, 'LineWidth', 2); legend('Estimated Velocity','Ground Truth','Location','NorthEast'); title('Velocity Along z Axis [m/s]'); xlabel('t [s]'); ylabel('Velocity [m/s]'); grid on; grid minor;xlim([67.9376, 70.0537]);
subplot(313), plot(time(1:end-N), x(9,:), 'LineWidth', 2); hold on; plot(time, AA1, 'LineWidth', 2); legend('Estimated Acceleration','Ground Truth','Location','SouthEast'); title('Acceleration Along z Axis [m/s^2]'); xlabel('t [s]'); ylabel('Acceleration [m/s^2]'); grid on; grid minor;xlim([67.9376, 70.0537]);

rmsErrorZ = rms(pos-PP1(2+N:end))
%
% figure
% subplot(9,1,1), plot(time, x(1,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along x axis'); grid on;
% subplot(9,1,2), plot(time, x(2,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along y axis'); grid on;
% subplot(9,1,3), plot(time, x(3,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(9,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along z axis'); grid on;
% subplot(9,1,4), plot(time, x(4,:), 'b', 'LineWidth', 2); %hold on; plot(time,, phi_gt(1,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along x axis'); grid on;
% subplot(9,1,5), plot(time, x(5,:), 'b', 'LineWidth', 2); %hold on; plot(time, phi_gt(2,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along x axis'); grid on;
% subplot(9,1,6), plot(time, x(6,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along y axis'); grid on;
% subplot(9,1,7), plot(time, x(7,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(4,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along y axis'); grid on;
% subplot(9,1,8), plot(time, x(8,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(5,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along z axis'); grid on;
% subplot(9,1,9), plot(time, x(9,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along z axis'); grid on;
