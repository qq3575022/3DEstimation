clc, clear, close all
data=readtable('2.csv','Delimiter', ',');  g=9.7953; load('RFxyz.mat');
%
figure
subplot(411), plot(t1, mag1)
subplot(412), plot(t2, mag2)
subplot(413), plot(t3, mag3)
subplot(414), plot(t4, mag4)
%%
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
xS  = find(abs(acc_time-107.99)<0.002);  xS = xS(1);    xE = find(abs(acc_time-111.984)<0.001);   xE = xE(1);
xxS = find(abs(gyro_time-107.99)<0.002); xxS = xxS(1); xxE = find(abs(gyro_time-111.984)<0.001); xxE = xxE(1);
xSS = find(abs(mag_time-107.99)<0.012);  xSS = xSS(1); xEE = find(abs(mag_time-111.984)<0.009);  xEE = xEE(1);
%time = unique(sort([accT; gyroT; magT]),'rows');

% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x = data_x(find_acc_data);
acc_data_y = -data_y(find_acc_data);
acc_data_z = data_z(find_acc_data)-g;

% AXY = [1, 0.0623, 0.0055; 0, 1, -0.0041; 0, 0, 1]*[0.9944, 0, 0; 0, 0.9999, 0; 0, 0, 0.9880]*([acc_data_x';acc_data_z';-acc_data_y'] + [0.1739; 0.0071; -0.2999]);
% 
% acc_data_x = AXY(1,:); acc_data_y = -AXY(3,:); acc_data_z = AXY(2,:);
% 

%%
gyro_data_x = data_x(find_gyro_data); gyro_data_y = data_y(find_gyro_data); gyro_data_z = data_z(find_gyro_data);
mag_data_x  = data_x(find_mag_data);   mag_data_y = data_y(find_mag_data);   mag_data_z = data_z(find_mag_data);

accx = acc_data_x(xS:xE) - mean(acc_data_x(xS:xE)) - 0.002*(acc_time(xS:xE) - acc_time(xS));
accy = acc_data_y(xS:xE) - mean(acc_data_y(xS:xE)) - 0.011*(acc_time(xS:xE) - acc_time(xS));
accz = acc_data_z(xS:xE) - mean(acc_data_z(xS:xE));
accT = acc_time(xS:xE);

gyrox = gyro_data_x(xxS:xxE) - mean(gyro_data_x(xxS:xxE));gyroy = gyro_data_y(xxS:xxE) - mean(gyro_data_y(xxS:xxE));gyroz = gyro_data_z(xxS:xxE) - mean(gyro_data_z(xxS:xxE));
gyroT = gyro_time(xxS:xxE);

magx = mag_data_x(xSS:xEE) - mean(mag_data_x(xSS:xEE));magy = mag_data_y(xSS:xEE) - mean(mag_data_y(xSS:xEE));magz = mag_data_z(xSS:xEE) - mean(mag_data_z(xSS:xEE));
magT = mag_time(xSS:xEE);
% ========== Start End Time Acceleration ========== 

time = unique(sort([accT; gyroT; magT]),'rows');
%
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
coord3(3,:) = PP3+1.09;

% Parameters of tag
Gt = 14.62;    % tag's antenna gain
X = 0.85;      % polarization mismatch
M = 4;         % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

x1 = [0,    0,  0.865];%[2.6256, 0.0889,0.858];%[-0.2,-0.2, 0.4];
x2 = [2.29, 0,  1.27];%[ 0.9, 0.8, 0.2];
x3 = [2.29,2.52, 0.865];%[ 0.8,-0.3,-0.2];
x4 = [0, 2.52,  1.27];

% Parameters of reader
PT = 1;         % reader's transmitted power
GT = 14.62;     % reader's trasmitter antenna gain 9.5dBi
GR = 14.62;     % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.00012; 

% phase cconcatenation
% global l1; global l2; global l3; l1 = 0;l2 = 0;l3 = 0;k = 1; 

% Get RSS and phase from each reader observing the moving tag
z = NaN(3,length(coord3)-1); z_prev = NaN(3,length(coord3)-1);

z_prev(1,:) = coord3(1,1:end-1); z(1,:) = coord3(1,2:end);% x coordinate
z_prev(2,:) = coord3(2,1:end-1); z(2,:) = coord3(2,2:end);% y coordinate
z_prev(3,:) = coord3(3,1:end-1); z(3,:) = coord3(3,2:end);% z coordinate

zdiff1 = z(1,:)-z_prev(1,:);
zdiff2 = z(2,:)-z_prev(2,:);
zdiff3 = z(3,:)-z_prev(3,:);

%
H1 = NaN(1,length(PP1));     phi1 = NaN(1,length(PP1));       phi_mu1 = NaN(1,length(PP1));
H2 = NaN(1,length(PP1));     phi2 = NaN(1,length(PP1));       phi_mu2 = NaN(1,length(PP1));
H3 = NaN(1,length(PP1));     phi3 = NaN(1,length(PP1));       phi_mu3 = NaN(1,length(PP1));
H4 = NaN(1,length(PP1));     phi4 = NaN(1,length(PP1));       phi_mu4 = NaN(1,length(PP1));

r_sim1 = NaN(1,length(PP1)); r_sim1_ = NaN(1,length(PP1)); rdot_sim1 = NaN(1,length(PP1)); rdot_sim1_ = NaN(1,length(PP1));  diff1 = NaN(1,length(PP1));
r_sim2 = NaN(1,length(PP1)); r_sim2_ = NaN(1,length(PP1)); rdot_sim2 = NaN(1,length(PP1)); rdot_sim2_ = NaN(1,length(PP1));  diff2 = NaN(1,length(PP1));
r_sim3 = NaN(1,length(PP1)); r_sim3_ = NaN(1,length(PP1)); rdot_sim3 = NaN(1,length(PP1)); rdot_sim3_ = NaN(1,length(PP1));  diff3 = NaN(1,length(PP1));
r_sim4 = NaN(1,length(PP1)); r_sim4_ = NaN(1,length(PP1)); rdot_sim4 = NaN(1,length(PP1)); rdot_sim4_ = NaN(1,length(PP1)); diff4 = NaN(1,length(PP1));

for k = 1:1:length(time)-1  
[H1(k+1),phi_mu1(k+1),r_sim1(k+1),r_sim1_(k+1),rdot_sim1(k+1),rdot_sim1_(k+1),diff1(k+1)] = noisysim(x1,f1,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,time(k+1)-time(k));
[H2(k+1),phi_mu2(k+1),r_sim2(k+1),r_sim2_(k+1),rdot_sim2(k+1),rdot_sim2_(k+1),diff2(k+1)] = noisysim(x2,f2,Gt,M,X,PT,GT,GR,R,sigma,2,k,z,z_prev,time(k+1)-time(k));
[H3(k+1),phi_mu3(k+1),r_sim3(k+1),r_sim3_(k+1),rdot_sim3(k+1),rdot_sim3_(k+1),diff3(k+1)] = noisysim(x3,f3,Gt,M,X,PT,GT,GR,R,sigma,3,k,z,z_prev,time(k+1)-time(k));   
[H4(k+1),phi_mu4(k+1),r_sim4(k+1),r_sim4_(k+1),rdot_sim4(k+1),rdot_sim4_(k+1),diff4(k+1)] = noisysim(x4,f4,Gt,M,X,PT,GT,GR,R,sigma,4,k,z,z_prev,time(k+1)-time(k));  
end

figure
subplot(4,1,1),plot(time, H1,'LineWidth',2);title('Simulated H from Reader $\#1$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
subplot(4,1,2),plot(time, H2,'LineWidth',2);title('Simulated H from Reader $\#2$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
subplot(4,1,3),plot(time, H3,'LineWidth',2);title('Simulated H from Reader $\#3$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')
subplot(4,1,4),plot(time, H4,'LineWidth',2);title('Simulated H from Reader $\#4$ in 2D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]')

figure
subplot(4,1,1),plot(t1+13.5167, mag1), hold on, plot(time, 0.01*H1,'LineWidth',2);
subplot(4,1,2),plot(t2+13.5167, mag2), hold on, plot(time, 0.01*H2,'LineWidth',2);
subplot(4,1,3),plot(t3+13.5167, mag3), hold on, plot(time, 0.01*H3,'LineWidth',2);
subplot(4,1,4),plot(t4+13.5167, mag4), hold on, plot(time, 0.01*H4,'LineWidth',2);
%%
%r_sim1_ = r_sim1+0.1*rand(1, length(r_sim1))-0.5*0.1; r_sim2_ = r_sim2+rand(1, length(r_sim2)); r_sim3_ = r_sim3+rand(1, length(r_sim3)); r_sim4_ = r_sim4+rand(1, length(r_sim4));
r_sim1(1) = r_sim1(2); r_sim2(1) = r_sim2(2); r_sim3(1) = r_sim3(2); r_sim4(1) = r_sim4(2); 
rdot_sim1(1) = rdot_sim1(2); rdot_sim2(1) = rdot_sim2(2); rdot_sim3(1) = rdot_sim3(2); rdot_sim4(1) = rdot_sim4(2); 
figure
subplot(8,1,1),plot(time, r_sim1_,'LineWidth',2);hold on;plot(time, r_sim1,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_1$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,2),plot(time, rdot_sim1_,'-','LineWidth',2);hold on;plot(time, rdot_sim1,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_1$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,3),plot(time, r_sim2_,'LineWidth',2);hold on;plot(time, r_sim2,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_2$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,4),plot(time, rdot_sim2_,'-','LineWidth',2);hold on;plot(time, rdot_sim2,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_2$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,5),plot(time, r_sim3_,'LineWidth',2);hold on;plot(time, r_sim3,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_3$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,6),plot(time, rdot_sim3_,'-','LineWidth',2);hold on;plot(time, rdot_sim3,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_3$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
subplot(8,1,7),plot(time, r_sim4_,'LineWidth',2);hold on;plot(time, r_sim4,'LineWidth',3);legend('Simulated','Ground truth');title('2D Radial Distance $R_4$','interpreter','latex');ylabel('radius [m]');xlabel('t [s]')
subplot(8,1,8),plot(time, rdot_sim4_,'-','LineWidth',2);hold on;plot(time, rdot_sim4,'-','LineWidth',2);legend('Simulated','Ground truth');title('2D Radial Velocity $\dot R_4$','interpreter','latex');ylabel('radial velocity [m/s]');xlabel('t [s]')
% Nonlinear - Extended kalman filter
% Load ground truth orientation, angular velocity, acc_x, acc_y, acc_z
%[phi,phi_gt,AXY,XX,Tmag,tdmag,Tacc,tdacc,T3,td3] = groundtruth3D();

mag1 = medfilt1(mag1,20);
mag2 = medfilt1(mag2,20);
mag3 = medfilt1(mag3,20);
mag4 = medfilt1(mag4,20);

time1 = t1 +time(1) - t1(1); 
time2 = t2 + time(1) - t2(1); 
time3 = t3+time(1) - t3(1);
mag11 = mag1; mag22 = mag2; mag33 = mag3; mag44 = mag4; 

%
phase11 = phase1; phase22 = phase2; phase33 = phase3; phase44 = phase4; 


phaseDiff1 = diff(phase11); phaseDiff2 = diff(phase22);  phaseDiff3 = diff(phase33);  phaseDiff4 = diff(phase44); 

phaseDiff1 = -medfilt1(phaseDiff1,200)*1/(42e-03);
phaseDiff2 = -medfilt1(phaseDiff2,200)*1/(42e-03);
phaseDiff3 = -medfilt1(phaseDiff3,200)*1/(42e-03);
phaseDiff4 = -medfilt1(phaseDiff4,200)*1/(42e-03);

%%
% 
% measr1 = 0.006*mag11.^(-0.5)+1.0; measr2 = 0.0005*mag22.^(-0.5)+1.79; 
% measr3 = 0.0042*mag33.^(-0.5)+0.9; measr4 = 0.0005*mag44.^(-0.5)+1.54;

measr1 = 0.006*mag11.^(-0.5)+1.0; measr2 = 0.0005*mag22.^(-0.5)+1.78; 
measr3 = 0.0042*mag33.^(-0.5)+0.98; measr4 = 0.0005*mag44.^(-0.5)+1.54;

% measr1 = 0.007*mag11.^(-0.5)+0.9; measr2 = 0.0005*mag22.^(-0.5)+1.78; 
% measr3 = 0.0042*mag33.^(-0.5)+0.94; measr4 = 0.0007*mag44.^(-0.5)+1.5;

% measr1 = medfilt1(measr1,10);
% measr2 = medfilt1(measr2,10);
% measr3 = medfilt1(measr3,10);
% measr4 = medfilt1(measr4,10);

measrdot1 = phaseDiff1; measrdot2 = phaseDiff2; measrdot3 = phaseDiff3; measrdot4 = phaseDiff4; 

dmeasr1 = zeros(1, length(time));dmeasr2 = zeros(1, length(time));
dmeasr3 = zeros(1, length(time));dmeasr4 = zeros(1, length(time));

dmeasrdot1 = zeros(1, length(time));  dmeasrdot2 = zeros(1, length(time)); 
dmeasrdot3 = zeros(1, length(time));  dmeasrdot4 = zeros(1, length(time)); 

for i = 3:1:length(time)
    difT = time(i) - time(1);
    index = min(round(difT*(4.8469*1000)), length(measrdot1));
    index2 = min(round(difT*(9.8607*1000)), length(measrdot3));
    dmeasr1(i) = measr1(index);
    dmeasr2(i) = measr2(index);
    dmeasr3(i) = measr3(index2);
    dmeasr4(i) = measr4(index);
    
    dmeasrdot1(i) = measrdot1(index);
    dmeasrdot2(i) = measrdot2(index);
    dmeasrdot3(i) = measrdot3(index2);
    dmeasrdot4(i) = measrdot4(index);
end

dmeasr1(1) = measr1(1);dmeasr2(1) = measr2(1);dmeasr3(1) = measr3(1);dmeasr4(1) = measr4(1);
dmeasrdot1(1) = measrdot1(1);dmeasrdot2(1) = measrdot2(1);dmeasrdot3(1) = measrdot3(1);dmeasrdot4(1) = measrdot4(1);
dmeasr1(2) = measr1(2);dmeasr2(2) = measr2(2);dmeasr3(2) = measr3(2);dmeasr4(2) = measr4(2);
dmeasrdot1(2) = measrdot1(2);dmeasrdot2(2) = measrdot2(2);dmeasrdot3(2) = measrdot3(2);dmeasrdot4(2) = measrdot4(2);


figure
subplot(411), plot(time,  dmeasr1, 'LineWidth', 1), hold on, plot(time, r_sim1, 'LineWidth', 1), title('Radial Distance $r1$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);
subplot(412), plot(time,  dmeasr2, 'LineWidth', 1), hold on, plot(time, r_sim2, 'LineWidth', 1),title('Radial Distance $r2$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);
subplot(413), plot(time,  dmeasr3, 'LineWidth', 1), hold on, plot(time, r_sim3, 'LineWidth', 1),title('Radial Distance $r3$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);
subplot(414), plot(time,  dmeasr4, 'LineWidth', 1), hold on, plot(time, r_sim4, 'LineWidth', 1),title('Radial Distance $r4$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);


figure
subplot(411), plot(time, dmeasrdot1); hold on; plot(time, rdot_sim1, 'LineWidth', 1), title('Radial Velocity $\dot r_1$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);
subplot(412), plot(time, dmeasrdot2); hold on; plot(time, rdot_sim2, 'LineWidth', 1), title('Radial Velocity $\dot r_2$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);
subplot(413), plot(time, dmeasrdot3); hold on; plot(time, rdot_sim3, 'LineWidth', 1), title('Radial Velocity $\dot r_3$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth','Location','SouthEast'), xlim([107.99, 111.984]);
subplot(414), plot(time, dmeasrdot4); hold on; plot(time, rdot_sim4, 'LineWidth', 1), title('Radial Velocity $\dot r_4$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);

%%
N = 3; i = 1; j = 1; k = 1;
% Initiate x, e and temproary variable
x = [0.5;0.5;0.5; 0.5;0.5;0.5; 1.1;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]*ones(1,length(time)-N);e = NaN(1,length(time)-N);
    
for m = 1:1:length(time)-N
    
    [y, i, j, k, index, len] = gety(dmeasr1, dmeasrdot1, dmeasr2, dmeasrdot2, dmeasr3, dmeasrdot3, dmeasr4, dmeasrdot4, magx, gyrox, magy, gyroy, magz, gyroz, accx, accy, accz, accT, gyroT, magT, time, i, j, k, m, N);

    x(:,m) = lsqnonlin(@(xx)getNLE(y, xx, N, index, len, time, m),[0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]);

    F = getF(time, m, m+N);
    x(:,m) = F*x(:,m);

end


%
[pos1, vel1] = doubleInte(diff(time(1:end-N)), x(3,2:end));
[pos2, vel2] = doubleInte(diff(time(1:end-N)), x(6,2:end));
[pos3, vel3] = doubleInte(diff(time(1:end-N)), x(9,2:end));

%%
figure
subplot(311), plot(time(1:end-N-1), 0.7*(pos1+1.03)+0.3*x(1,2:end), 'LineWidth', 2); hold on; plot(time, PP1+1.03, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','NorthEast'); title('Position Along x Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([107.99, 111.984])
subplot(312), plot(time(1:end-N-1), 0.7*vel1+0.3*x(2,2:end), 'LineWidth', 2); hold on; plot(time, VV1, 'LineWidth', 2); legend('Estimated Velocity','Ground Truth'); title('Velocity Along x Axis [m/s]'); xlabel('t [s]'); ylabel('Velocity [m/s]'); grid on; grid minor;xlim([107.99, 111.984])%xlim([3.5, 7.5])
subplot(313), plot(time(1:end-N), x(3,:), 'LineWidth', 2); hold on; plot(time, AA1, 'LineWidth', 2); legend('Estimated Acceleration','Ground Truth'); title('Acceleration Along x Axis [m/s^2]'); xlabel('t [s]'); ylabel('Acceleration [m/s^2]'); grid on; grid minor;xlim([107.99, 111.984])


figure
subplot(311), plot(time(1:end-N-1), 0.7*(pos2+1.31)+0.3*x(4,2:end), 'LineWidth', 2); hold on; plot(time, PP2+1.31, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','NorthEast'); title('Position Along y Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([107.99, 111.984])
subplot(312), plot(time(1:end-N-1), 0.7*vel2+0.3*x(5,2:end), 'LineWidth', 2); hold on; plot(time, VV2, 'LineWidth', 2); legend('Estimated Velocity','Ground Truth'); title('Velocity Along y Axis [m/s]'); xlabel('t [s]'); ylabel('Velocity [m/s]'); grid on; grid minor;xlim([107.99, 111.984])%xlim([3.5, 7.5])
subplot(313), plot(time(1:end-N), x(6,:), 'LineWidth', 2); hold on; plot(time, AA2, 'LineWidth', 2); legend('Estimated Acceleration','Ground Truth'); title('Acceleration Along y Axis [m/s^2]'); xlabel('t [s]'); ylabel('Acceleration [m/s^2]'); grid on; grid minor;xlim([107.99, 111.984])

figure
subplot(311), plot(time(1:end-N-1), 0.7*(pos3+1.03)+0.3*x(7,2:end), 'LineWidth', 2); hold on; plot(time, PP3+1.03, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','SouthEast'); title('Position Along z Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([107.99, 111.984])
subplot(312), plot(time(1:end-N-1), 0.7*vel3+0.3*x(8,2:end), 'LineWidth', 2); hold on; plot(time, VV3, 'LineWidth', 2); legend('Estimated Velocity','Ground Truth'); title('Velocity Along z Axis [m/s]'); xlabel('t [s]'); ylabel('Velocity [m/s]'); grid on; grid minor;xlim([107.99, 111.984])%xlim([3.5, 7.5])
subplot(313), plot(time(1:end-N), x(9,:), 'LineWidth', 2); hold on; plot(time, AA3, 'LineWidth', 2); legend('Estimated Acceleration','Ground Truth'); title('Acceleration Along z Axis [m/s^2]'); xlabel('t [s]'); ylabel('Acceleration [m/s^2]'); grid on; grid minor;xlim([107.99, 111.984])


%
figure
subplot(311), plot(time(1:end-N-1), 0.7*(pos1+1.03)+0.3*x(1,2:end)); hold on; plot(time, PP1+1.03, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','NorthEast'); title('Position Along x Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([107.99, 111.984])
subplot(312), plot(time(1:end-N-1), 0.7*(pos2+1.31)+0.3*x(4,2:end)); hold on; plot(time, PP2+1.31, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','NorthEast'); title('Position Along y Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([107.99, 111.984])
subplot(313), plot(time(1:end-N-1), 0.7*(pos3+1.03)+0.3*x(7,2:end)); hold on; plot(time, PP2+1.03, 'LineWidth', 2); legend('Estimated Position','Ground Truth','Location','NorthEast'); title('Position Along z Axis [m]'); xlabel('t [s]'); ylabel('Position [m]'); grid on; grid minor;xlim([107.99, 111.984])

%%
rmsErrorX = rms(0.7*(pos1+1.03)+0.3*x(1,2:end)-PP1(1:end-N-1)-1.03)
rmsErrorY = rms(0.7*(pos2+1.31)+0.3*x(4,2:end)-PP2(1:end-N-1)-1.31)
rmsErrorZ = rms(0.7*(pos3+1.03)+0.3*x(7,2:end)-PP3(1:end-N-1)-1.03)
%%
%
figure
subplot(9,1,1), plot(time(1:end-N), x(1,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along x axis'); grid on;
subplot(9,1,2), plot(time(1:end-N), x(2,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along y axis'); grid on;
subplot(9,1,3), plot(time(1:end-N), x(3,:), 'b', 'LineWidth', 2); %hold on; plot(time, XX(9,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along z axis'); grid on;
subplot(9,1,4), plot(time(1:end-N), x(4,:), 'b', 'LineWidth', 2); %hold on; plot(time,, phi_gt(1,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along x axis'); grid on;
subplot(9,1,5), plot(time(1:end-N), x(5,:), 'b', 'LineWidth', 2); %hold on; plot(time, phi_gt(2,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along x axis'); grid on;
subplot(9,1,6), plot(time(1:end-N), x(6,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(3,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along y axis'); grid on;
subplot(9,1,7), plot(time(1:end-N), x(7,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(4,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along y axis'); grid on;
subplot(9,1,8), plot(time(1:end-N), x(8,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(5,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along z axis'); grid on;
subplot(9,1,9), plot(time(1:end-N), x(9,:), 'b', 'LineWidth', 2); %hold on; plot(td3, phi_gt(6,:), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along z axis'); grid on;
