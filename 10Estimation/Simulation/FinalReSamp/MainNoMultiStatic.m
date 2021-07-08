clc, clear, close all

%
% IMU 
% data=readtable('2.csv','Delimiter', ',');  g=9.7953;
% 
% % Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
% sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
% sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
% sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);
% 
% imu_time   = table2array(data(:,1));               data_x     = table2array(data(:,3)); 
% acc_time   = imu_time(find_acc_data);              data_y     = table2array(data(:,4));    
% acc_time   = (acc_time - acc_time(1))/1000000000;  data_z     = table2array(data(:,5));
% 
% gyro_time   = imu_time(find_gyro_data);
% gyro_time   = (gyro_time - gyro_time(1))/1000000000;
% 
% xS = find(abs(acc_time-107.99)<0.002); xS = xS(1);
% xE = find(abs(acc_time-111.984)<0.002); xE = xE(1);
% 
% yS = find(abs(acc_time-107.99)<0.002); yS = yS(1);
% yE = find(abs(acc_time-111.984)<0.002); yE = yE(1);
% 
% zS = find(abs(acc_time-107.99)<0.002); zS = zS(1);
% zE = find(abs(acc_time-111.984)<0.002); zE = zE(1);
% 
% 
% accT = acc_time(xS:xE); gyroT = gyro_time(xS:xE);
% %
% % ==========   Raw Acceleration Data  ========== 
% % ----acc_data_x----
% % ----acc_data_y----
% % ----acc_data_z----
% acc_data_x = data_x(find_acc_data);
% acc_data_y = data_y(find_acc_data);
% acc_data_z = data_z(find_acc_data);
% 
% 
% figure
% subplot(311),plot(acc_time, acc_data_x,'LineWidth',3)%, xlim([0,40]);,%hold on, plot(acc_time, A3_sim(1,:),'LineWidth',2),  hold on, plot(acc_time, A3(1,:),'LineWidth',2),  xlim([0,10]); xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along x Axe Based on Error Model')
% subplot(312),plot(acc_time, acc_data_z,'LineWidth',3)%, xlim([0,40]);,%hold on, plot(acc_time, -A3_sim(2,:),'LineWidth',2), hold on, plot(acc_time, -A3(2,:),'LineWidth',2), xlim([45,55]); xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along y Axe Based on Error Model')
% subplot(313),plot(acc_time, acc_data_y,'LineWidth',3)%, xlim([0,40]);,%hold on,plot(acc_time, A3_sim(3,:),'LineWidth',2),hold on,plot(acc_time, A3(3,:),'LineWidth',2),xlim([90,100]); xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along z Axe Based on Error Model')
% %
% %
% tdx = acc_time(xS:xE) - acc_time(xS);
% tdy = acc_time(yS:yE) - acc_time(yS);
% tdz = acc_time(zS:zE) - acc_time(zS);
% %
% [PP1, VV1, AA1] = groundtruth1Dx2(tdx);
% [PP2, VV2, AA2] = groundtruth1Dy2(tdy);
% [PP3, VV3, AA3] = groundtruth1Dz2(tdz);
% 
% % Start from 0
% PP1 = PP1 - PP1(1);
% PP2 = PP2 - PP2(1);
% PP3 = PP3 - PP3(1);
% 
% % 1D - 3D coordinates
% % td3
% 
% A3 =  zeros(3, length(acc_time));
% % x
% A3(1,xS+1:length(AA1) + xS) = AA1;
% %y
% A3(2,yS+1:length(AA2) + yS) = AA2;
% %z
% A3(3,zS+1:length(AA3) + zS) = AA3;
% 
% % A3(3,length(AA3)*3+1:length(AA3)*4) = flip(AA3);
% 
% %     Ta = [1 0.0936 0.0621;
% %           0 1      0.0329;
% %           0 0      1];
% %     Ka =  diag([1.001, 0.9812, 0.9441]);
% %     
% %     TKa = pinv(Ta*Ka);
% %     
% %     errorA = TKa*A3;
% %     errorA3 = errorA - [0.0961;0.0928;0.0947];
% 
% 
% % ========== Filtered Acceleration Data ========== 
% % ----  accx  ----
% % ----  accy  ----
% % ----  accz  ----
% accx = medfilt1(acc_data_x,50);
% accy = medfilt1(acc_data_y,50) - g;
% accz = medfilt1(acc_data_z,50);
% 
% 
% figure
% subplot(311),plot(acc_time, accx,'LineWidth',3)%, xlim([0,40]);,%hold on, plot(acc_time, A3_sim(1,:),'LineWidth',2),  hold on, plot(acc_time, A3(1,:),'LineWidth',2),  xlim([0,10]); xlabel('time [s]'), ylabel('x axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along x Axe Based on Error Model')
% subplot(312),plot(acc_time, accz,'LineWidth',3)%, xlim([0,40]);,%hold on, plot(acc_time, -A3_sim(2,:),'LineWidth',2), hold on, plot(acc_time, -A3(2,:),'LineWidth',2), xlim([45,55]); xlabel('time [s]'), ylabel('y axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along y Axe Based on Error Model')
% subplot(313),plot(acc_time, accy,'LineWidth',3)%, xlim([0,40]);,%hold on,plot(acc_time, A3_sim(3,:),'LineWidth',2),hold on,plot(acc_time, A3(3,:),'LineWidth',2),xlim([90,100]); xlabel('time [s]'), ylabel('z axis ($m/s^2$)','interpreter','latex'),legend('Measurement', 'Simulated', 'Ground Truth with Bias');title('Simulated Accelerations Along z Axe Based on Error Model')
% 
% %
% 
% A3 = A3 + [-0.6373;-0.25;0.05]*ones(1, length(acc_time));
% A3_sim = A3 + [0;-0.000;0]*acc_time';
% 
% A3_sim = awgn(A3_sim,30,'measured');
% 
% acc = NaN(3, length(accT));
% acc(1,:) = A3_sim(1,xS:xE); acc(2,:) = A3_sim(2,xS:xE); acc(3,:) = A3_sim(3,xS:xE);
% 
% angle = zeros(3, length(gyroT)) + 0.01*randn(1, length(gyroT)); 
% gyro = zeros(3, length(gyroT))+ 0.01*randn(1, length(gyroT)); 
% 
% figure
% subplot(311), plot(gyroT, angle(1,:));
% subplot(312), plot(gyroT, angle(2,:));
% subplot(313), plot(gyroT, angle(3,:));
% 
% figure
% subplot(311), plot(gyroT, gyro(1,:));
% subplot(312), plot(gyroT, gyro(2,:));
% subplot(313), plot(gyroT, gyro(3,:));
% 
% figure
% subplot(311), plot(accT, acc(1,:));
% subplot(312), plot(accT, acc(2,:));
% subplot(313), plot(accT, acc(3,:));
% save('SimIMU.mat', 'acc', 'accT', 'angle', 'gyro', 'gyroT');

%%
% =================================================================== Load Data ======================================================================
% ----------- Position of Four Readers ---------
x1 = [0,    0,    0.865];  
x2 = [2.29, 0,    1.27];   
x3 = [2.29, 2.52, 0.865]; 
x4 = [0,    2.52, 1.27];

% -------------- Time and Coordinates ----------
% time:     IMU Measurement
% coord3:   Ground Truth of 3D Coordinates
% z,z_prev: 3D coordinates for Simulation
[magD12, magD22, magD32, magD42, time] = getMeas();% Measurement Magnitude of Length 130854
%%

% -------------- Time and Coordinates ----------
% time:     IMU Measurement
% coord3:   Ground Truth of 3D Coordinates
% z,z_prev: 3D coordinates for Simulation

[coord3, radial, z, z_prev, H1, H2, H3, H4, H1_, H2_, H3_, H4_, r_sim, r_sim1, r_sim2, r_sim3, r_sim4, r_sim1_, r_sim2_, r_sim3_, r_sim4_, r_meas, r_meas1, r_meas2, r_meas3, r_meas4, r_phase, rphase1, rphase2, rphase3, rphase4, rdot_sim1, rdot_sim2, rdot_sim3, rdot_sim4, rdot_sim1_, rdot_sim2_, rdot_sim3_, rdot_sim4_, phase1, phase2, phase3, phase4, phigt1_1, phigt1_2, phigt1_3, phigt1_4, phigt2_1, phigt2_2, phigt2_3, phigt2_4, phi1_1, phi1_2, phi1_3, phi1_4, phi2_1, phi2_2, phi2_3, phi2_4, phi3_1, phi3_2, phi3_3, phi3_4, phi4_1, phi4_2, phi4_3, phi4_4] = get3Dcoord(x1, x2, x3, x4, time);
%%
% Start and End Index of 3D Motion
yST  = find(abs(time-107.99)<0.002);   yST = yST(1)-1;
yET  = find(abs(time-111.984)<0.002);  yET = yET(1)+1;
%%
% =================================================================== Simulation ======================================================================

% Parameters of reader
Gt = 0.6*1.462*3/4;    % tag's antenna gain
X  = 0.85;             % polarization mismatch
M  = sqrt(4);          % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

% Parameters of reader
PT = 1;         % reader's transmitted power
GT = sqrt(1.462);     % reader's trasmitter antenna gain 9.5dBi
GR = sqrt(1.462);     % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.0012; 

offset11 = 0; offset21 = 0; offset31 = 0; offset41 = 0;offset12 = 0; offset22 = 0; offset32 = 0; offset42 = 0;
offset13 = 0; offset23 = 0; offset33 = 0; offset43 = 0;offset14 = 0; offset24 = 0; offset34 = 0; offset44 = 0;
offset51 = 0; offset52 = 0; offset53 = 0; offset54 = 0;offset61 = 0; offset62 = 0; offset63 = 0; offset64 = 0; 

rng(2,'twister');

for k = yST-1:1:yET  

[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1),r_meas1(k+1),rphase1(k+1),rdot_sim1(k+1),rdot_sim1_(k+1),phigt1_1(k+1),phigt2_1(k+1),phi1_1(k+1),phi2_1(k+1),offset11,offset21,offset31,offset41,offset51,offset61] = noisysimNoMultiStatic(x1,f1,Gt,M,X,PT,0.7*GT,GR,R,sigma,k,z,z_prev,phigt1_1(k),phi1_1(k),time(k)-time(k-1),magD12(k),offset11,offset21,offset31,offset41,offset51,offset61);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1),r_meas2(k+1),rphase2(k+1),rdot_sim2(k+1),rdot_sim2_(k+1),phigt1_2(k+1),phigt2_2(k+1),phi1_2(k+1),phi2_2(k+1),offset12,offset22,offset32,offset42,offset52,offset62] = noisysimNoMultiStatic(x2,f2,Gt,M,X,PT,  7*GT,GR,R,sigma,k,z,z_prev,phigt1_2(k),phi1_2(k),time(k)-time(k-1),magD22(k),offset12,offset22,offset32,offset42,offset52,offset62);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1),r_meas3(k+1),rphase3(k+1),rdot_sim3(k+1),rdot_sim3_(k+1),phigt1_3(k+1),phigt2_3(k+1),phi1_3(k+1),phi2_3(k+1),offset13,offset23,offset33,offset43,offset53,offset63] = noisysimNoMultiStatic(x3,f3,Gt,M,X,PT,    GT,GR,R,sigma,k,z,z_prev,phigt1_3(k),phi1_3(k),time(k)-time(k-1),magD32(k),offset13,offset23,offset33,offset43,offset53,offset63);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1),r_meas4(k+1),rphase4(k+1),rdot_sim4(k+1),rdot_sim4_(k+1),phigt1_4(k+1),phigt2_4(k+1),phi1_4(k+1),phi2_4(k+1),offset14,offset24,offset34,offset44,offset54,offset64] = noisysimNoMultiStatic(x4,f4,Gt,M,X,PT,0.5*GT,GR,R,sigma,k,z,z_prev,phigt1_4(k),phi1_4(k),time(k)-time(k-1),magD42(k),offset14,offset24,offset34,offset44,offset54,offset64);  

rphase1(k+1) = rphase1(k+1)+0;
rphase2(k+1) = rphase2(k+1)+0;
rphase3(k+1) = rphase3(k+1)+0;
rphase4(k+1) = rphase4(k+1)+0;

[r_sim(1,k+1),r_sim(2,k+1),r_sim(3,k+1)]   = trilateration(x1, x2, x3, x4, r_sim1(k+1), r_sim2(k+1), r_sim3(k+1),  r_sim4(k+1));
[r_phase(1,k+1),r_phase(2,k+1),r_phase(3,k+1)] = trilateration(x1, x2, x3, x4, rphase1(k+1),rphase2(k+1),rphase3(k+1),rphase4(k+1));
[r_meas(1,k+1),r_meas(2,k+1),r_meas(3,k+1)]  = trilateration(x1, x2, x3, x4, r_meas1(k+1),r_meas2(k+1),r_meas3(k+1),r_meas4(k+1));

end

rtime = time(yST:yET); 

r1 = rphase1(yST:yET); r2 = rphase2(yST:yET); r3 = rphase3(yST:yET); r4 = rphase4(yST:yET);
rdot1 = rdot_sim1(yST:yET); rdot2 = rdot_sim2(yST:yET); rdot3 = rdot_sim3(yST:yET); rdot4 = rdot_sim4(yST:yET);

save('SimRFPhase1.mat', 'r1', 'r2', 'r3', 'r4', 'rdot1', 'rdot2', 'rdot3', 'rdot4','rtime');
r1 = r_sim1(yST:yET); r2 = r_sim2(yST:yET); r3 = r_sim3(yST:yET); r4 = r_sim4(yST:yET);

save('SimRFMag1.mat', 'r1', 'r2', 'r3', 'r4', 'rdot1', 'rdot2', 'rdot3', 'rdot4','rtime');

%%

%% ==================================================================== Results ======================================================================
% -------------- Radial Distances RMS: Reader1 Reader2 Reader3 Reader4 ----------
fprintf('================================== radial distance RMS ================================== ')
rmeas  = 1000*[rms(r_meas1(yST:yET) - radial(1,yST:yET)), rms(r_meas2(yST:yET) - radial(2,yST:yET)), rms(r_meas3(yST:yET) - radial(3,yST:yET)), rms(r_meas4(yST:yET) - radial(4,yST:yET))]
rmag   = 1000*[rms(r_sim1(yST:yET)  - radial(1,yST:yET)), rms(r_sim2(yST:yET)  - radial(2,yST:yET)), rms(r_sim3(yST:yET) - radial(3,yST:yET)),  rms(r_sim4(yST:yET)  - radial(4,yST:yET))]
rphase = 1000*[rms(rphase1(yST:yET) - radial(1,yST:yET)), rms(rphase2(yST:yET) - radial(2,yST:yET)), rms(rphase3(yST:yET) - radial(3,yST:yET)), rms(rphase4(yST:yET) - radial(4,yST:yET))]

fprintf('================================== x y z triangulation RMS in total ================================== ')
meas  = 1000*sqrt(rms(coord3(1,yST:yET) - r_meas(1,yST:yET)).^2  + rms(coord3(2,yST:yET) - r_meas(2,yST:yET)).^2  + rms(coord3(3,yST:yET) - r_meas(3,yST:yET)).^2)
mag   = 1000*sqrt(rms(coord3(1,yST:yET) - r_sim(1,yST:yET)).^2   + rms(coord3(2,yST:yET) - r_sim(2,yST:yET)).^2   + rms(coord3(3,yST:yET) - r_sim(3,yST:yET)).^2)
phase = 1000*sqrt(rms(coord3(1,yST:yET) - r_phase(1,yST:yET)).^2 + rms(coord3(2,yST:yET) - r_phase(2,yST:yET)).^2 + rms(coord3(3,yST:yET) - r_phase(3,yST:yET)).^2)
%%
fprintf('================================== x y z triangulation RMS separately ================================== ')
meas  = 1000*[rms(coord3(1,yST:yET) - r_meas(1,yST:yET)), rms(coord3(2,yST:yET) - r_meas(2,yST:yET)),  rms(coord3(3,yST:yET) - r_meas(3,yST:yET))]
mag   = 1000*[rms(coord3(1,yST:yET) - r_sim(1,yST:yET)),  rms(coord3(2,yST:yET) - r_sim(2,yST:yET)),   rms(coord3(3,yST:yET) - r_sim(3,yST:yET))]
phase = 1000*[rms(coord3(1,yST:yET) - r_phase(1,yST:yET)),rms(coord3(2,yST:yET) - r_phase(2,yST:yET)), rms(coord3(3,yST:yET) - r_phase(3,yST:yET))]

%%
fprintf('================================== radial distance mean value ================================== ')
rmeas  = 100*[mean(r_meas1(yST:yET) - radial(1,yST:yET)), mean(r_meas2(yST:yET) - radial(2,yST:yET)), mean(r_meas3(yST:yET) - radial(3,yST:yET)), mean(r_meas4(yST:yET) - radial(4,yST:yET))]
rmag   = 100*[mean(r_sim1(yST:yET)  - radial(1,yST:yET)), mean(r_sim2(yST:yET)  - radial(2,yST:yET)), mean(r_sim3(yST:yET) - radial(3,yST:yET)),  mean(r_sim4(yST:yET)  - radial(4,yST:yET))]
rphase = 100*[mean(rphase1(yST:yET) - radial(1,yST:yET)), mean(rphase2(yST:yET) - radial(2,yST:yET)), mean(rphase3(yST:yET) - radial(3,yST:yET)), mean(rphase4(yST:yET) - radial(4,yST:yET))]

fprintf('================================== x y z triangulation Mean Error in total ================================== ')
meas  = 100*sqrt(mean(coord3(1,yST:yET) - r_meas(1,yST:yET)).^2  + mean(coord3(2,yST:yET) - r_meas(2,yST:yET)).^2  + mean(coord3(3,yST:yET) - r_meas(3,yST:yET)).^2)
mag   = 100*sqrt(mean(coord3(1,yST:yET) - r_sim(1,yST:yET)).^2   + mean(coord3(2,yST:yET) - r_sim(2,yST:yET)).^2   + mean(coord3(3,yST:yET) - r_sim(3,yST:yET)).^2)
phase = 100*sqrt(mean(coord3(1,yST:yET) - r_phase(1,yST:yET)).^2 + mean(coord3(2,yST:yET) - r_phase(2,yST:yET)).^2 + mean(coord3(3,yST:yET) - r_phase(3,yST:yET)).^2)

%%
% -------------- Measurement Radial Distances: Reader1 Reader2 Reader3 Reader4 ----------

figure
subplot(4,1,1),plot(time(yST:yET), r_meas1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), r_meas2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), r_meas3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), r_meas4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

% -------------- Magnitude Radial Distances: Reader1 Reader2 Reader3 Reader4 ----------
%%
figure
subplot(4,1,1),plot(time(yST:yET), r_sim1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([1, 3]);legend('Simulation','Ground truth');
subplot(4,1,2),plot(time(yST:yET), r_sim2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([0.8, 2.8]);legend('Simulation','Ground truth');
subplot(4,1,3),plot(time(yST:yET), r_sim3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([0.5, 2.5]);legend('Simulation','Ground truth');
subplot(4,1,4),plot(time(yST:yET), r_sim4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([0.8, 2.8]);legend('Simulation','Ground truth');

% -------------- Two Frequencies Radial Distances: Reader1 Reader2 Reader3 Reader4 ----------
%%
figure
subplot(4,1,1),plot(time(yST:yET), rphase1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,2),plot(time(yST:yET), rphase2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,3),plot(time(yST:yET), rphase3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,4),plot(time(yST:yET), rphase4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');

%%
figure
subplot(3,1,1),plot(time(yST:yET), r_phase(1,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(1,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(3,1,2),plot(time(yST:yET), r_phase(2,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(2,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(3,1,3),plot(time(yST:yET), r_phase(3,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(3,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');

%%
figure
subplot(3,1,1),plot(time(yST:yET), r_sim(1,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(1,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(3,1,2),plot(time(yST:yET), r_sim(2,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(2,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(3,1,3),plot(time(yST:yET), r_sim(3,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(3,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
%%
figure
subplot(4,1,1),plot(time(yST:yET), rdot_sim1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), rdot_sim1_(yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,2),plot(time(yST:yET), rdot_sim2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), rdot_sim2_(yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,3),plot(time(yST:yET), rdot_sim3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), rdot_sim3_(yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,4),plot(time(yST:yET), rdot_sim4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), rdot_sim4_(yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
