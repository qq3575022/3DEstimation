clc, clear, close all

x1 = [0,    0,    0.865];  
x2 = [2.29, 0,    1.27];   
x3 = [2.29, 2.52, 0.865]; 
x4 = [0,    2.52, 1.27];

load('RF.mat'); 
load('IMU.mat')
%% -------------- Time and Coordinates ----------
% time:     IMU Measurement
% coord3:   Ground Truth of 3D Coordinates
% z,z_prev: 3D coordinates for Simulation

itime = unique(sort([accT; gyroT]),'rows');
%
time = unique(sort([rtime; itime]),'rows');

[coord3, radial] = get3Dcoord(x1, x2, x3, x4, time);
%%
N = 5; i = 1; j = 1; k = 1; l = 1;
% Initiate x, e and temproary variable
x = [0.5;0.5;0.5; 0.5;0.5;0.5; 1.1;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]*ones(1,length(time)-N);e = NaN(1,length(time)-N);
factor = 3;

for m = 1:1:length(time)-N
    
    [y, i, j, l, index, len] = gety(r1, rdot1, r2, rdot2, r3+0.05, rdot3, r4, rdot4, angle(1,:), gyro(1,:), angle(2,:), gyro(2,:), angle(3,:), gyro(3,:), acc(1,:), acc(2,:), acc(3,:), accT, gyroT, rtime, time, i, j, l, m, N, factor);

    x(:,m) = lsqnonlin(@(xx)getNLE(y, xx, N, index, len, time, m, factor),[0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]);

    F = getF(time, m, m+N);
    x(:,m) = F*x(:,m);

end


%
[pos1, vel1] = doubleInte(diff(time(1:end-N)), x(3,2:end));
[pos2, vel2] = doubleInte(diff(time(1:end-N)), x(6,2:end));
[pos3, vel3] = doubleInte(diff(time(1:end-N)), x(9,2:end));

%%
figure
subplot(9,1,1), plot(time(1:end-N), x(1,:),  'LineWidth', 3); hold on; plot(time, coord3(1,1:end), 'r', 'LineWidth', 3); legend('estimated', 'ground truth','location','southeast'); title('Estimated position along x axis'); grid on; grid minor; ylabel('pos  [m]');xlim([108, 112])
subplot(9,1,2), plot(time(1:end-N), x(2,:),  'LineWidth', 3); hold on; plot(time, coord3(4,1:end), 'r', 'LineWidth', 3); legend('estimated', 'ground truth'); title('Estimated velocity along x axis'); grid on;grid minor; ylabel('vel  [m/s]');xlim([108, 112])
subplot(9,1,3), plot(time(1:end-N), x(3,:),  'LineWidth', 3); hold on; plot(time, coord3(7,1:end), 'r', 'LineWidth', 3); legend('estimated', 'ground truth'); title('Estimated acceleration along x axis'); grid on;grid minor; ylabel('acc  [m/s^2]');xlim([108, 112])
subplot(9,1,4), plot(time(1:end-N), x(4,:),  'LineWidth', 3); hold on; plot(time, coord3(2,1:end), 'r', 'LineWidth', 3); legend('estimated', 'ground truth','location','southeast'); title('Estimated position along y axis'); grid on; grid minor; ylabel('pos  [m]');xlim([108, 112])
subplot(9,1,5), plot(time(1:end-N), x(5,:),  'LineWidth', 3); hold on; plot(time, coord3(5,1:end), 'r', 'LineWidth', 3); legend('estimated', 'ground truth'); title('Estimated velocity along y axis'); grid on;grid minor; ylabel('vel  [m/s]');xlim([108, 112])
subplot(9,1,6), plot(time(1:end-N), x(6,:),  'LineWidth', 3); hold on; plot(time, coord3(8,1:end), 'r', 'LineWidth', 3); legend('estimated', 'ground truth'); title('Estimated acceleration along y axis'); grid on;grid minor; ylabel('acc  [m/s^2]');xlim([108, 112])
subplot(9,1,7), plot(time(1:end-N), x(7,:),  'LineWidth', 3); hold on; plot(time, coord3(3,1:end), 'r', 'LineWidth', 3); legend('estimated', 'ground truth','location','southeast'); title('Estimated position along z axis'); grid on; grid minor; ylabel('pos  [m]');xlim([108, 112])
subplot(9,1,8), plot(time(1:end-N), x(8,:),  'LineWidth', 3); hold on; plot(time, coord3(6,1:end), 'r', 'LineWidth', 3); legend('estimated', 'ground truth'); title('Estimated velocity along z axis'); grid on;grid minor; ylabel('vel  [m/s]');xlim([108, 112])
subplot(9,1,9), plot(time(1:end-N), x(9,:),  'LineWidth', 3); hold on; plot(time, coord3(9,1:end), 'r', 'LineWidth', 3); legend('estimated', 'ground truth'); title('Estimated acceleration along z axis'); grid on;grid minor; ylabel('acc  [m/s^2]');xlim([108, 112])
xlabel('t [s]')
%%
rmsErrorX = 1000*rms(x(1,:)-coord3(1,1:end-N));
rmsErrorY = 1000*rms(x(4,:)-coord3(2,1:end-N));
rmsErrorZ = 1000*rms(x(7,:)-coord3(3,1:end-N));
rms3D = sqrt(rmsErrorX^2 + rmsErrorY^2 + rmsErrorZ^2)
rms2D = sqrt(rmsErrorX^2 + rmsErrorY^2)
rmsErrorZ