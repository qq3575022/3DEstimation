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
% yS = find(abs(acc_time-46.13)<0.002); yS = yS(1);yE = find(abs(acc_time-48.5)<0.002); yE = yE(1);
% yyS = find(abs(gyro_time-46.13)<0.002);  yyS = yyS(1);yyE = find(abs(gyro_time-48.5)<0.002); yyE = yyE(1);
% ySS = find(abs(mag_time-46.13)<0.009);   ySS = ySS(1);yEE = find(abs(mag_time-48.5)<0.009); yEE = yEE(1);

yS = find(abs(acc_time-57.9476)<0.002); yS = yS(1);
yE = find(abs(acc_time-61.0160)<0.002); yE = yE(1);

yyS = find(abs(gyro_time-57.9476)<0.002);  yyS = yyS(1);
yyE = find(abs(gyro_time-61.0160)<0.002); yyE = yyE(1);

ySS = find(abs(mag_time-57.9476)<0.012);   ySS = ySS(1);
yEE = find(abs(mag_time-61.0160)<0.012); yEE = yEE(1);

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

accx = acc_data_x(yS:yE) - mean(acc_data_x(yS:yE));accy = acc_data_y(yS:yE) - mean(acc_data_y(yS:yE));accz = acc_data_z(yS:yE) - mean(acc_data_z(yS:yE)) + 0.008;
accT = acc_time(yS:yE);

gyrox = gyro_data_x(yyS:yyE) - mean(gyro_data_x(yyS:yyE));gyroy = gyro_data_y(yyS:yyE) - mean(gyro_data_y(yyS:yyE));gyroz = gyro_data_z(yyS:yyE) - mean(gyro_data_z(yyS:yyE));

angle(1,:) = angle(1,:) - mean(angle(1,:)); 
angle(2,:) = angle(2,:) - mean(angle(2,:)); 
angle(3,:) = angle(3,:) - mean(angle(3,:)); 


gyroT = gyro_time(yyS:yyE);


magx = mag_data_x(ySS:yEE) - mean(mag_data_x(ySS:yEE));magy = mag_data_y(ySS:yEE) - mean(mag_data_y(ySS:yEE));magz = mag_data_z(ySS:yEE) - mean(mag_data_z(ySS:yEE));
magT = mag_time(ySS:yEE);
% ========== Start End Time Acceleration ========== 

time = unique(sort([accT; gyroT; magT]),'rows');

[PP1, VV1, AA1] = groundtruth1Dy2(time - time(1));

% Start from 0
PP1 = PP1 - PP1(1);
%
%time = unique(sort([accT; gyroT; magT]),'rows');
%
%[PP1, VV1, AA1] = groundtruth1Dx(time - time(1)); PP1 = PP1 - PP1(1);

N = 5; i = 1; j = 1; k = 1;
% Initiate x, e and temproary variable
x = [0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]*ones(1,length(time)-N);e = NaN(1,length(time)-N);
    
for m = 1:1:length(time)-N
    
    [y, i, j, k, index, len] = gety(angle(1,:), gyrox, angle(2,:), gyroy, angle(3,:), gyroz, accx, accy, accz, accT, gyroT, magT, time, i, j, k, m);

    x(:,m) = lsqnonlin(@(xx)getNLE(y, xx, N, index, len, time, m),[0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]);

    F = getF(time, m, m+1);
    x(:,m) = F*x(:,m);

    
end
%
[pos, vel] = doubleInte(diff(time(1:end-N)), x(6,2:end));

%%
figure
subplot(311), plot(time(2:end-N), pos, 'LineWidth', 3); hold on; plot(time, PP1, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','NorthEast'); title('Position Along y Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([57.9176, 61.0160]);
subplot(312), plot(time(2:end-N), vel, 'LineWidth', 2); hold on; plot(time, VV1, 'LineWidth', 2); legend('Estimated Velocity','Ground Truth','Location','SouthEast'); title('Velocity Along y Axis [m/s]'); xlabel('t [s]'); ylabel('Velocity [m/s]'); grid on; grid minor;xlim([57.9176, 61.0160]);
subplot(313), plot(time(1:end-N), x(6,:), 'LineWidth', 2); hold on; plot(time, AA1, 'LineWidth', 2); legend('Estimated Acceleration','Ground Truth','Location','SouthEast'); title('Acceleration Along y Axis [m/s^2]'); xlabel('t [s]'); ylabel('Acceleration [m/s^2]'); grid on; grid minor;xlim([57.9176, 61.0160]);

%%
rmsErrorY = rms(pos-PP1(N+2:end))
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