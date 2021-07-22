clc, clear, close all

x1 = [0,    0,    0.865];  
x2 = [2.29, 0,    1.27];   
x3 = [2.29, 2.52, 0.865]; 
x4 = [0,    2.52, 1.27];

load('RF.mat'); 
%load('IMU.mat')
% -------------- Time and Coordinates ----------
% time:     IMU Measurement
% coord3:   Ground Truth of 3D Coordinates
% z,z_prev: 3D coordinates for Simulation


[coord3, radial] = get3Dcoord(x1, x2, x3, x4, rtime);
%%
N = 3; i = 1; j = 1; k = 1; l = 1;
% Initiate x, e and temproary variable
x = [0.5;0.5;0.5; 0.5;0.5;0.5; 1.1;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]*ones(1,length(rtime)-N);e = NaN(1,length(rtime)-N);
factor = 1;

for m = 1:1:length(rtime)-N
    
    [y,  l, index, len] = gety(r1, r2, r3, r4, rtime, l, m, N, factor);

    x(:,m) = lsqnonlin(@(xx)getNLE(y, xx, N, index, len, rtime, m, factor),[0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]);

    F = getF(rtime, m, m+N);
    x(:,m) = F*x(:,m);

end


%
[pos1, vel1] = doubleInte(diff(rtime(1:end-N)), x(3,2:end));
[pos2, vel2] = doubleInte(diff(rtime(1:end-N)), x(6,2:end));
[pos3, vel3] = doubleInte(diff(rtime(1:end-N)), x(9,2:end));

%%
figure
subplot(9,1,1), plot(rtime(1:end-N), x(1,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(1,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along x axis'); grid on;
subplot(9,1,2), plot(rtime(1:end-N), x(2,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(4,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along y axis'); grid on;
subplot(9,1,3), plot(rtime(1:end-N), x(3,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(7,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along z axis'); grid on;
subplot(9,1,4), plot(rtime(1:end-N), x(4,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(2,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along x axis'); grid on;
subplot(9,1,5), plot(rtime(1:end-N), x(5,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(5,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along x axis'); grid on;
subplot(9,1,6), plot(rtime(1:end-N), x(6,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(8,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along y axis'); grid on;
subplot(9,1,7), plot(rtime(1:end-N), x(7,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(3,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along y axis'); grid on;
subplot(9,1,8), plot(rtime(1:end-N), x(8,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(6,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along z axis'); grid on;
subplot(9,1,9), plot(rtime(1:end-N), x(9,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(9,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along z axis'); grid on;


%%
rmsErrorX = 1000*rms(x(1,:)-coord3(1,1:end-N))
rmsErrorY = 1000*rms(x(4,:)-coord3(2,1:end-N))
rmsErrorZ = 1000*rms(x(7,:)-coord3(3,1:end-N))