clc, clear, close all
data=readtable('z1.csv','Delimiter', ',');  g=9.7953;

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

time   = table2array(data(:,1));  data_x = table2array(data(:,3));  data_y = table2array(data(:,4));  data_z = table2array(data(:,5));

acc_time   = time(find_acc_data);  acc_time   = (acc_time - acc_time(1))/1000000000;                
gyro_time  = time(find_gyro_data); gyro_time  = (gyro_time - gyro_time(1))/1000000000;
mag_time   = time(find_mag_data);  mag_time   = (mag_time - mag_time(1))/1000000000;

time  = unique(sort([acc_time; gyro_time; mag_time]),'rows');

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x  = data_x(find_acc_data);   acc_data_y = data_y(find_acc_data)-g; acc_data_z = data_z(find_acc_data);
gyro_data_x = data_x(find_gyro_data);  gyro_data_y = data_y(find_gyro_data); gyro_data_z = data_z(find_gyro_data);
mag_data_x  = data_x(find_mag_data);   mag_data_y = data_y(find_mag_data);   mag_data_z = data_z(find_mag_data);

%%
% 
acc_data_x = medfilt1(acc_data_x,80);
acc_data_y = medfilt1(acc_data_y,80);
acc_data_z = medfilt1(acc_data_z,80);
% 

figure

subplot(911), plot(acc_time(5682/2+800:end), acc_data_x(5682/2+800:end),'LineWidth',2),grid on 
grid minor
ylabel('Acceleration [m/s^2]'); title('Acceleration x');


subplot(912), plot(acc_time(5682/2+800:end), acc_data_z(5682/2+800:end),'LineWidth',2),grid on 
grid minor;ylabel('Acceleration [m/s^2]');title('Acceleration y');

subplot(913), plot(acc_time(5682/2+800:end), acc_data_y(5682/2+800:end),'LineWidth',2),grid on 
grid minor;ylabel('Acceleration [m/s^2]');title('Acceleration z');


subplot(914), plot(gyro_time, gyro_data_x,'LineWidth',2),grid on 
grid minor
hold on, plot(gyro_time, zeros(1, length(gyro_data_x)),'LineWidth',2);ylabel('Angular Velocity [rad/s]');title('Angular Velocity x [rad/s]');

subplot(915), plot(gyro_time, gyro_data_y,'LineWidth',2);grid on 
grid minor
hold on, plot(gyro_time, zeros(1, length(gyro_data_x)),'LineWidth',2);ylabel('Angular Velocity [rad/s]');title('Angular Velocity y [rad/s]');

subplot(916), plot(gyro_time, gyro_data_z,'LineWidth',2);grid on 
grid minor
hold on, plot(gyro_time, zeros(1, length(gyro_data_x)),'LineWidth',2);ylabel('Angular Velocity [rad/s]');title('Angular Velocity z [rad/s]');

subplot(917), plot(mag_time, mag_data_x,'LineWidth',2);grid on 
grid minor
ylabel('Magnitude Strength');title('Magnitude Strength');

subplot(918), plot(mag_time, mag_data_y,'LineWidth',2);grid on 
grid minor
ylabel('Magnitude Strength');title('Magnitude Strength');

subplot(919), plot(mag_time, mag_data_z,'LineWidth',2);grid on 
grid minor
ylabel('Magnitude Strength');title('Magnitude Strength');
