clc, clear, close all
data=readtable('2.csv','Delimiter', ',');  g=9.7953;

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

acc_time   = table2array(data(:,1));               data_x     = table2array(data(:,3)); 
acc_time   = acc_time(find_acc_data);              data_y     = table2array(data(:,4));    
acc_time   = (acc_time - acc_time(1))/1000000000;  data_z     = table2array(data(:,5));

xS = find(abs(acc_time-47.8957)<0.002); xS = xS(1);
xE = find(abs(acc_time-51.452)<0.002); xE = xE(1);

yS = find(abs(acc_time-57.9176)<0.002); yS = yS(1);
yE = find(abs(acc_time-61.0160)<0.002); yE = yE(1);

zS = find(abs(acc_time-67.9376)<0.002); zS = zS(1);
zE = find(abs(acc_time-70.0537)<0.002); zE = zE(1);

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x = data_x(find_acc_data);
acc_data_y = data_y(find_acc_data);
acc_data_z = data_z(find_acc_data);


figure
subplot(311),plot(acc_time, acc_data_x,'LineWidth',3)%, xlim([0,40]);,%hold on, plot(acc_time, A3_sim(1,:),'LineWidth',2),  hold on, plot(acc_time, A3(1,:),'LineWidth',2),  xlim([0,10]); xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along x Axe Based on Error Model')
subplot(312),plot(acc_time, acc_data_z,'LineWidth',3)%, xlim([0,40]);,%hold on, plot(acc_time, -A3_sim(2,:),'LineWidth',2), hold on, plot(acc_time, -A3(2,:),'LineWidth',2), xlim([45,55]); xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(313),plot(acc_time, acc_data_y,'LineWidth',3)%, xlim([0,40]);,%hold on,plot(acc_time, A3_sim(3,:),'LineWidth',2),hold on,plot(acc_time, A3(3,:),'LineWidth',2),xlim([90,100]); xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along z Axe Based on Error Model')
%
%
tdx = acc_time(xS:xE) - acc_time(xS);
tdy = acc_time(yS:yE) - acc_time(yS);
tdz = acc_time(zS:zE) - acc_time(zS);
%
[PP1, VV1, AA1] = groundtruth1Dx2(tdx);
[PP2, VV2, AA2] = groundtruth1Dy2(tdy);
[PP3, VV3, AA3] = groundtruth1Dz2(tdz);

% Start from 0
PP1 = PP1 - PP1(1);
PP2 = PP2 - PP2(1);
PP3 = PP3 - PP3(1);

% 1D - 3D coordinates
% td3

A3 =  zeros(3, length(acc_time));
% x
A3(1,xS+1:length(AA1) + xS) = AA1;
%y
A3(2,yS+1:length(AA2) + yS) = AA2;
%z
A3(3,zS+1:length(AA3) + zS) = AA3;

% A3(3,length(AA3)*3+1:length(AA3)*4) = flip(AA3);

%     Ta = [1 0.0936 0.0621;
%           0 1      0.0329;
%           0 0      1];
%     Ka =  diag([1.001, 0.9812, 0.9441]);
%     
%     TKa = pinv(Ta*Ka);
%     
%     errorA = TKa*A3;
%     errorA3 = errorA - [0.0961;0.0928;0.0947];


% ========== Filtered Acceleration Data ========== 
% ----  accx  ----
% ----  accy  ----
% ----  accz  ----
accx = medfilt1(acc_data_x,50);
accy = medfilt1(acc_data_y,50) - g;
accz = medfilt1(acc_data_z,50);


figure
subplot(311),plot(acc_time, accx,'LineWidth',3)%, xlim([0,40]);,%hold on, plot(acc_time, A3_sim(1,:),'LineWidth',2),  hold on, plot(acc_time, A3(1,:),'LineWidth',2),  xlim([0,10]); xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along x Axe Based on Error Model')
subplot(312),plot(acc_time, accz,'LineWidth',3)%, xlim([0,40]);,%hold on, plot(acc_time, -A3_sim(2,:),'LineWidth',2), hold on, plot(acc_time, -A3(2,:),'LineWidth',2), xlim([45,55]); xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(313),plot(acc_time, accy,'LineWidth',3)%, xlim([0,40]);,%hold on,plot(acc_time, A3_sim(3,:),'LineWidth',2),hold on,plot(acc_time, A3(3,:),'LineWidth',2),xlim([90,100]); xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along z Axe Based on Error Model')

A3 = A3 + [-0.6373;-0.25;0.05]*ones(1, length(acc_time));
A3_sim = A3 + [0;-0.0000;0]*acc_time';

A3_sim = awgn(A3_sim,30,'measured');

[mean(A3_sim(1,:) - A3(1,:)), rms(A3_sim(1,:) - A3(1,:)), sqrt(rms(A3_sim(1,:) - A3(1,:)))]
[mean(A3_sim(2,:) - A3(2,:)), rms(A3_sim(2,:) - A3(2,:)), sqrt(rms(A3_sim(2,:) - A3(2,:)))]
[mean(A3_sim(3,:) - A3(3,:)), rms(A3_sim(3,:) - A3(3,:)), sqrt(rms(A3_sim(3,:) - A3(3,:)))]

[mean(accx' - A3(1,:)), rms(accx' - A3(1,:)), sqrt(rms(accx' - A3(1,:)))]
[mean(accz' + A3(2,:)), rms(accz' + A3(2,:)), sqrt(rms(accz' + A3(2,:)))]
[mean(accy' - A3(3,:)), rms(accy' - A3(3,:)), sqrt(rms(accy' - A3(3,:)))]

figure
subplot(611),plot(acc_time, A3_sim(1,:),'LineWidth',2), xlim([47.9457, 51.452]);%hold on, plot(acc_time, accx,'LineWidth',2), xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along x Axe Based on Error Model')
subplot(612),plot(acc_time, accx,'LineWidth',2),xlim([47.9457, 51.452]), %ylim([-0.5, 0]);%hold on, plot(acc_time, accz,'LineWidth',2), xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(613),plot(acc_time, A3_sim(2,:),'LineWidth',2),xlim([57.8076, 61.0161]); %hold on, plot(acc_time, accz,'LineWidth',2), xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(614),plot(acc_time, accz,'LineWidth',2),xlim([57.8076, 61.0161]), %ylim([-0.6, 0.1]); %hold on, plot(acc_time, accz,'LineWidth',2), xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(615),plot(acc_time, A3_sim(3,:),'LineWidth',2),xlim([67.9376, 70.0537]); %hold on, plot(acc_time, accy,'LineWidth',2), xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along z Axe Based on Error Model')
subplot(616),plot(acc_time(10:length(acc_time)), accy(10:length(acc_time)),'LineWidth',2), ylim([-0.2, 0.2]),xlim([67.9376, 70.0537]); %hold on, plot(acc_time, accy,'LineWidth',2), xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),title('Simulated Accelerations Along z Axe Based on Error Model')

figure
subplot(311),plot(acc_time, accx,'LineWidth',3),hold on, plot(acc_time, A3_sim(1,:),'LineWidth',2),  hold on, plot(acc_time, A3(1,:),'LineWidth',2),  xlim([47, 52]); xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along x Axe Based on Error Model')
subplot(312),plot(acc_time, accz,'LineWidth',3),hold on, plot(acc_time, -A3_sim(2,:),'LineWidth',2), hold on, plot(acc_time, -A3(2,:),'LineWidth',2), xlim([56, 63]); xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along y Axe Based on Error Model')
subplot(313),plot(acc_time(10:length(acc_time)), accy(10:length(acc_time)),'LineWidth',3),hold on,plot(acc_time, A3_sim(3,:),'LineWidth',2),hold on,plot(acc_time, A3(3,:),'LineWidth',2),xlim([66, 72]); xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along z Axe Based on Error Model')

