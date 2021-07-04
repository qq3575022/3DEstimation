function [coord3, radial] = get3Dcoord(x1, x2, x3, x4, time)



%xyz
xyzSSS = find(abs(time-107.99)<0.009);   xyzSSS = xyzSSS(1); xyzEEE = find(abs(time-111.984)<0.009);      xyzEEE = xyzEEE(1);


[PPxyz1, VVxyz1, AAxyz1] = groundtruth1Dx2(time-time(1)); % move along x y z simulanteously
[PPxyz2, VVxyz2, AAxyz2] = groundtruth1Dy2(time-time(1)); % move along x y z simulanteously
[PPxyz3, VVxyz3, AAxyz3] = groundtruth1Dz2(time-time(1)); % move along x y z simulanteously

%
% 3D coordinates
coord3 = zeros(9, length(time));
%
% x
coord3(1,:) = PPxyz1 + 1.03;
%coord3(1,xyzSSS+1+length(PPxyz1):end) = PPxyz1(end);
%
%y
coord3(2,:) = PPxyz2 + 1.31;
%coord3(2,xyzSSS+1+length(PPxyz2):end) = PPxyz1(end);
%z
coord3(3,:) = PPxyz3 + 1.03;
%coord3(3,xyzSSS+1+length(PPxyz3):end) = PPxyz3(end);
%===============================================================
%=====================================
% x
coord3(4,:) = VVxyz1;
%coord3(4,xyzSSS+1+length(PPxyz1):end) = VVxyz1(end);
%

%y
coord3(5,:) = VVxyz2;
%coord3(5,xyzSSS+1+length(PPxyz2):end) = VVxyz2(end);

%z
coord3(6,:) = VVxyz3;
%coord3(6,xyzSSS+1+length(PPxyz3):end) = VVxyz3(end);

%===============================================================
%=====================================
% x
coord3(7,:) = AAxyz1;
%coord3(7,xyzSSS+1+length(PPxyz1):end) = AAxyz1(end);
%
%y
coord3(8,:) = AAxyz2;
%coord3(8,xyzSSS+1+length(PPxyz2):end) = AAxyz2(end);

%z
coord3(9,:) = AAxyz3;
%coord3(9,xyzSSS+1+length(PPxyz3):end) = AAxyz3(end);
%%
% get z H1 r_sim
% Get RSS and phase from each reader observing the moving tag

radial = NaN(4, length(coord3));

for i = 1:1:length(coord3)
    radial(1,i) = sqrt((coord3(1,i) - x1(1))^2+(coord3(2,i) - x1(2))^2+(coord3(3,i) - x1(3))^2);
    radial(2,i) = sqrt((coord3(1,i) - x2(1))^2+(coord3(2,i) - x2(2))^2+(coord3(3,i) - x2(3))^2);
    radial(3,i) = sqrt((coord3(1,i) - x3(1))^2+(coord3(2,i) - x3(2))^2+(coord3(3,i) - x3(3))^2);
    radial(4,i) = sqrt((coord3(1,i) - x4(1))^2+(coord3(2,i) - x4(2))^2+(coord3(3,i) - x4(3))^2);
end

end


% % ==========   Raw Acceleration Data  ========== 
% % ----acc_data_x----
% % ----acc_data_y----
% % ----acc_data_z----
% acc_data_x  = data_x(find_acc_data);   acc_data_y = data_y(find_acc_data)-g; acc_data_z = data_z(find_acc_data);
% gyro_data_x = data_x(find_gyro_data);  gyro_data_y = data_y(find_gyro_data); gyro_data_z = data_z(find_gyro_data);
% mag_data_x  = data_x(find_mag_data);   mag_data_y = data_y(find_mag_data);   mag_data_z = data_z(find_mag_data);
% 
% % 
% acc_data_x = medfilt1(acc_data_x,80);
% acc_data_y = medfilt1(acc_data_y,80);
% acc_data_z = medfilt1(acc_data_z,80);
% % 
% %x
% accxX = acc_data_x(xS:xE) - mean(acc_data_x(xS:xE)); accyX = acc_data_y(xS:xE) - mean(acc_data_y(xS:xE)); acczX = acc_data_z(xS:xE) - mean(acc_data_z(xS:xE));
% accTX = acc_time(xS:xE);
% 
% gyroxX = gyro_data_x(xxS:xxE) - mean(gyro_data_x(xxS:xxE)); gyroyX = gyro_data_y(xxS:xxE) - mean(gyro_data_y(xxS:xxE)); gyrozX = gyro_data_z(xxS:xxE) - mean(gyro_data_z(xxS:xxE));
% gyroTX = gyro_time(xxS:xxE);
% 
% magxX = mag_data_x(xSS:xEE) - mean(mag_data_x(xSS:xEE)); magyX = mag_data_y(xSS:xEE) - mean(mag_data_y(xSS:xEE)); magzX = mag_data_z(xSS:xEE) - mean(mag_data_z(xSS:xEE));
% magTX = mag_time(xSS:xEE);
% % timeX
% timeX = unique(sort([accTX; gyroTX; magTX]),'rows');
% 
% %y
% accxY = acc_data_x(yS:yE) - mean(acc_data_x(yS:yE)); accyY = acc_data_y(yS:yE) - mean(acc_data_y(yS:yE)); acczY = acc_data_z(yS:yE) - mean(acc_data_z(yS:yE));
% accTY = acc_time(yS:yE);
% 
% gyroxY = gyro_data_x(yyS:yyE) - mean(gyro_data_x(yyS:yyE)); gyroyY = gyro_data_y(yyS:yyE) - mean(gyro_data_y(yyS:yyE)); gyrozY = gyro_data_z(yyS:yyE) - mean(gyro_data_z(yyS:yyE));
% gyroTY = gyro_time(yyS:yyE);
% 
% magxY = mag_data_x(ySS:yEE) - mean(mag_data_x(ySS:yEE)); magyY = mag_data_y(ySS:yEE) - mean(mag_data_y(ySS:yEE)); magzY = mag_data_z(ySS:yEE) - mean(mag_data_z(ySS:yEE));
% magTY = mag_time(ySS:yEE);
% % timeY
% timeY = unique(sort([accTY; gyroTY; magTY]),'rows');
% 
% %z
% accxZ = acc_data_x(zS:zE) - mean(acc_data_x(zS:zE)); accyZ = acc_data_y(zS:zE) - mean(acc_data_y(zS:zE)); acczZ = acc_data_z(zS:zE) - mean(acc_data_z(zS:zE));
% accTZ = acc_time(zS:zE);
% 
% gyroxZ = gyro_data_x(zzS:zzE) - mean(gyro_data_x(zzS:zzE)); gyroyZ = gyro_data_y(zzS:zzE) - mean(gyro_data_y(zzS:zzE)); gyrozZ = gyro_data_z(zzS:zzE) - mean(gyro_data_z(zzS:zzE));
% gyroTZ = gyro_time(zzS:zzE);
% 
% magxZ = mag_data_x(zSS:zEE) - mean(mag_data_x(zSS:zEE)); magyZ = mag_data_y(zSS:zEE) - mean(mag_data_y(zSS:zEE)); magzZ = mag_data_z(zSS:zEE) - mean(mag_data_z(zSS:zEE));
% magTZ = mag_time(zSS:zEE);
% % timeZ
% timeZ = unique(sort([accTZ; gyroTZ; magTZ]),'rows');
% 
% %xyz back
% accxbZ = acc_data_x(bS1:bE1) - mean(acc_data_x(bS1:bE1)); accybZ = acc_data_y(bS1:bE1) - mean(acc_data_y(bS1:bE1)); acczbZ = acc_data_z(bS1:bE1) - mean(acc_data_z(bS1:bE1));
% accTbZ = acc_time(bS1:bE1);
% 
% gyroxbZ = gyro_data_x(bbS1:bbE1) - mean(gyro_data_x(bbS1:bbE1)); gyroybZ = gyro_data_y(bbS1:bbE1) - mean(gyro_data_y(bbS1:bbE1)); gyrozbZ = gyro_data_z(bbS1:bbE1) - mean(gyro_data_z(bbS1:bbE1));
% gyroTbZ = gyro_time(bbS1:bbE1);
% 
% magxbZ = mag_data_x(bSS1:bEE1) - mean(mag_data_x(bSS1:bEE1)); magybZ = mag_data_y(bSS1:bEE1) - mean(mag_data_y(bSS1:bEE1)); magzbZ = mag_data_z(bSS1:bEE1) - mean(mag_data_z(bSS1:bEE1));
% magTbZ = mag_time(bSS1:bEE1);
% 
% % timeZ
% timebZ = unique(sort([accTbZ; gyroTbZ; magTbZ]),'rows');
% 
% % xyz
% accxZ = acc_data_x(xyzS:xyzE) - mean(acc_data_x(xyzS:xyzE)); accyZ = acc_data_y(xyzS:xyzE) - mean(acc_data_y(xyzS:xyzE)); acczZ = acc_data_z(xyzS:xyzE) - mean(acc_data_z(xyzS:xyzE));
% accTXYZ = acc_time(xyzS:xyzE);
% 
% gyroxZ = gyro_data_x(xyzzS:xyzzE) - mean(gyro_data_x(xyzzS:xyzzE)); gyroyZ = gyro_data_y(xyzzS:xyzzE) - mean(gyro_data_y(xyzzS:xyzzE)); gyrozZ = gyro_data_z(xyzzS:xyzzE) - mean(gyro_data_z(xyzzS:xyzzE));
% gyroTXYZ = gyro_time(xyzzS:xyzzE);
% 
% magxZ = mag_data_x(xyzSS:xyzEE) - mean(mag_data_x(xyzSS:xyzEE)); magyZ = mag_data_y(xyzSS:xyzEE) - mean(mag_data_y(xyzSS:xyzEE)); magzZ = mag_data_z(xyzSS:xyzEE) - mean(mag_data_z(xyzSS:xyzEE));
% magTXYZ = mag_time(xyzSS:xyzEE);
% 
% 
% % timeXYZ
% timeXYZ = unique(sort([accTXYZ; gyroTXYZ; magTXYZ]),'rows');
% 
% % Acceleration
% 
% % x then y finally z
% % [P1, V1, A1] = groundtruth1Dx(accTX - accTX(1));
% % [P2, V2, A2] = groundtruth1Dy(accTY - accTY(1));
% % [P3, V3, A3] = groundtruth1Dz(accTZ - accTZ(1));
% % 
% % % xyz together
% % [Pxyz1, Vxyz1, Axyz1] = groundtruth1Dx2(accTXYZ-accTXYZ(1));
% % [Pxyz2, Vxyz2, Axyz2] = groundtruth1Dy2(accTXYZ-accTXYZ(1));
% % [Pxyz3, Vxyz3, Axyz3] = groundtruth1Dz2(accTXYZ-accTXYZ(1));
% % 
% % %
% % AA =  zeros(3, length(acc_time));
% % 
% % % x
% % AA(1,xS+1:length(A1) + xS) = A1; AA(1,xyzS+1:length(Axyz1) + xyzS) = Axyz1; %AA(1,bS1+1:length(Axyz1) + bS1) = -Axyz1;%A3(1,xyzSback+1:length(AA1) + xyzSback) = -AA1; A3(1,xyzS+1:length(AA1) + xyzS) = AA1;
% % %y
% % AA(2,yS+1:length(A2) + yS) = A2; AA(2,xyzS+1:length(Axyz2) + xyzS) = Axyz2; %AA(2,bS1+1:length(Axyz2) + bS1) = -Axyz2;%A3(2,xyzSback+1:length(AA2) + xyzSback) = -AA2; A3(2,xyzS+1:length(AA2) + xyzS) = AA2;
% % %z
% % AA(3,zS+1:length(A3) + zS) = A3; AA(3,xyzS+1:length(Axyz3) + xyzS) = Axyz3; %AA(3,bS1+1:length(Axyz3) + bS1) = -Axyz3;%A3(3,xyzSback+1:length(AA3) + xyzSback) = -AA3; A3(3,xyzS+1:length(AA3) + xyzS) = AA3;
% % %time = unique(sort([accT; gyroT; magT]),'rows');
% 
% % All Measurement