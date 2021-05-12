%%
clc, clear, close all
%%
% Nonlinear - Extended kalman filter
% Load measurement orientation, angular velocity, acc_x, acc_y, acc_z

data=readtable('2.csv','Delimiter', ',');  g=9.7953; 

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

time   = table2array(data(:,1));               
data_x = table2array(data(:,3));  data_y = table2array(data(:,4));  data_z = table2array(data(:,5));

tdacc   = time(find_acc_data);  tdacc   = (tdacc - tdacc(1))/1000000000;                
tdgyro  = time(find_gyro_data); tdgyro  = (tdgyro - tdgyro(1))/1000000000;
tdmag   = time(find_mag_data);  tdmag   = (tdmag - tdmag(1))/1000000000;

td3 = unique(sort([tdacc; tdgyro; tdmag]),'rows'); T3 = diff(td3);

Tmag = diff(tdmag);
Tgyro = diff(tdgyro);
Tacc = diff(tdacc);

acc_data_x = data_x(find_acc_data);
acc_data_y = data_y(find_acc_data);
acc_data_z = -data_z(find_acc_data);

gyro_data_x = data_x(find_gyro_data); gyro_data_y = data_y(find_gyro_data); gyro_data_z = data_z(find_gyro_data);
mag_data_x  = data_x(find_mag_data);   mag_data_y = data_y(find_mag_data);   mag_data_z = data_z(find_mag_data);

mag = zeros(3, length(mag_data_x)); mag(1,:) = mag_data_x;  mag(2,:) = mag_data_y;  mag(3,:) = mag_data_z; 
gyro = zeros(3, length(gyro_data_x)); gyro(1,:) = gyro_data_x;  gyro(2,:) = gyro_data_y;  gyro(3,:) = gyro_data_z; 
AXY = zeros(3, length(acc_data_x)); AXY(1,:) = acc_data_x;  AXY(2,:) = acc_data_y;  AXY(3,:) = acc_data_z; 


% figure
% subplot(311), plot(tdacc, AXY(1,:))
% subplot(312), plot(tdacc, AXY(2,:))
% subplot(313), plot(tdacc, AXY(3,:))
% 
% figure
% subplot(311), plot(tdgyro, gyro(1,:))
% subplot(312), plot(tdgyro, gyro(2,:))
% subplot(313), plot(tdgyro, gyro(3,:))
% 
% figure
% subplot(311), plot(tdmag, mag(1,:))
% subplot(312), plot(tdmag, mag(2,:))
% subplot(313), plot(tdmag, mag(3,:))


%%
%[mag,gyro,AXY,Tmag,tdmag,Tgyro,tdgyro,Tacc,tdacc,T3,td3] = meas3D();

AXY = [1, 0.0623, 0.0055; 0, 1, -0.0041; 0, 0, 1]*[0.9944, 0, 0; 0, 0.9999, 0; 0, 0, 0.9880]*(AXY + [0.1739; 0.0071; -0.2999]);

%AXY(3,:) = AXY(3,:);% - 9.7932;

%%
% Measurement Vector [acc_x, acc_y, acc_z, orientation, angular velocity]

% State vector [q0 q1 q2 q3 B sigma]
x = [0.5, 0.5, 0.5, 0.5]'*ones(1,length(tdgyro));

% Intial P matrix
P = 0.8*eye(length(x(:,1)));
% covariance for w
Q = diag([0.3, 0.3, 0.3, 0.3]);

% covariance for v
RR = [var(AXY(1,:)), var(AXY(2,:)), var(AXY(3,:)), var(mag(1,:)), var(mag(2,:)), var(mag(3,:))]; 

angle = zeros(3,length(tdgyro));

n1 = 1;
n2 = 1;
n3 = 1;

for m = 1:1:length(tdgyro)

    %[index, i, j, k] = getLength(i, j, k, tdmag, tdgyro, tdacc);
 
    [y, index, n1, n2] = getyNPVA(tdgyro(m), tdacc, tdmag, mag, gyro, AXY, n1, n2);
    
    if index == 3
        
        index
    end
    
    if index == -1
        x(:,m) = x(:,m-1);
    else
        F = getF(gyro,Tgyro, m);

        if m == 1
            x(:,m) = x(:,m);
        else
            x(:,m) = F*x(:,m-1);
        end
        
        P = F*P*F' + Q;
        H = getJacoN(y,x(:,m),index);

        R = getR(RR, index);
        % back tracking
        G = P*H'/(H*P*H' + R);%lambda
        % e
        x(:,m) = x(:,m) + G*(y - getH(x(:,m), index));
        P = (eye(length(x(:,1))) - G*H)*P;
    end
    
    q0 = x(1,m)/sqrt(x(1,m)^2+x(2,m)^2+x(3,m)^2+x(4,m)^2); q1 = x(2,m)/sqrt(x(1,m)^2+x(2,m)^2+x(3,m)^2+x(4,m)^2); 
    q2 = x(3,m)/sqrt(x(1,m)^2+x(2,m)^2+x(3,m)^2+x(4,m)^2); q3 = x(4,m)/sqrt(x(1,m)^2+x(2,m)^2+x(3,m)^2+x(4,m)^2); 
    
%     angle(1,m) = -atan2(q2*q3-q0*q1, q0*q0+q3*q3-0.5);
%     angle(2,m) = asin(q1*q3+q0*q2);
%     angle(3,m) = -atan2(q1*q2-q0*q3, q0*q0+q1*q1-0.5);
    
    %angle(1,m) = atan2(2*(q2*q3 + q0*q1), q0*q0+q3*q3-q1*q1-q2*q2);
    angle(1,m) = atan2(2*(q2*q3 + q0*q1), q0*q0-q1*q1-q2*q2+q3*q3);%2*(q1*q1+q2*q2)-1);atan2(2*(q2*q3 + q0*q1), 2*(q1*q1+q2*q2)-1);
    %angle(2,m) = asin(-q1*q3+q0*q2);
    angle(2,m) = asin(2*(-q1*q3+q0*q2));
    %angle(3,m) = atan2(q1*q2 + q0*q3, q0*q0+q1*q1-q2*q2-q3*q3);
    angle(3,m) = atan2(2*(q1*q2 + q0*q3), 2*(q2*q2+q3*q3)-1);
    
%     if angle(1,m) > pi/2 || angle(1,m) < -pi/2
%         if angle(2,m) > 0.15
%             angle(2,m) = pi - angle(2,m);
%         elseif angle(2,m) < -0.25
%             angle(2,m) = -pi - angle(2,m);
%         end
%     end

end

% for m = 2:1:length(tdgyro)    
% %     if angle(1,m) - angle(1,m-1) > 5
% %         angle(1,m) = angle(1,m)-2*pi;
% %     elseif angle(1,m) - angle(1,m-1) < -5
% %         angle(1,m) = angle(1,m)+2*pi;
% %     end
%     
%     if abs(angle(2,m) - angle(2,m-1)) > 4.5
%         angle(2,m) = angle(2,m-1);
%     end
%     
%     if abs(angle(3,m) - angle(3,m-1)) > 4.5
%         angle(3,m) = angle(3,m-1);
%     end
% end

figure
subplot(3,1,1), plot(tdgyro, angle(1,:), 'LineWidth', 3); title('Estimated orientation along x axis -Roll $\theta_x$','interpreter','latex'); grid on; ylabel('$\theta_x$ [rad]','interpreter','latex')
subplot(3,1,2), plot(tdgyro, angle(2,:), 'LineWidth', 3); title('Estimated orientation along y axis -Pitch $\theta_y$','interpreter','latex'); grid on;ylabel('$\theta_y$ [rad]','interpreter','latex')
subplot(3,1,3), plot(tdgyro, angle(3,:), 'LineWidth', 3); title('Estimated orientation along z axis -Yaw $\theta_z$','interpreter','latex'); grid on; xlabel('time [s]');ylabel('$\theta_z$ [rad]','interpreter','latex')

%%
save('tdgyro.mat', 'tdgyro')
save('angle.mat', 'angle')