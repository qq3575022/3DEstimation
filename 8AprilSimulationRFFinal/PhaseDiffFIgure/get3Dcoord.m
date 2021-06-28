function [coord3, radial, z, z_prev, H1, H2, H3, H4, H1_, H2_, H3_, H4_, r_sim, r_sim1, r_sim2, r_sim3, r_sim4, r_sim1_, r_sim2_, r_sim3_, r_sim4_, r_meas, r_meas1, r_meas2, r_meas3, r_meas4, rphase, rphase1, rphase2, rphase3, rphase4, rdot_sim1, rdot_sim2, rdot_sim3, rdot_sim4, rdot_sim1_, rdot_sim2_, rdot_sim3_, rdot_sim4_, phase1, phase2, phase3, phase4, phigt1_1, phigt1_2, phigt1_3, phigt1_4, phigt2_1, phigt2_2, phigt2_3, phigt2_4, phi1_1, phi1_2, phi1_3, phi1_4, phi2_1, phi2_2, phi2_3, phi2_4, phi3_1, phi3_2, phi3_3, phi3_4, phi4_1, phi4_2, phi4_3, phi4_4] = get3Dcoord(x1, x2, x3, x4, time)

% D01 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader1_6.bin', 1000000000, 1);
% D1 = D01;%(85000:445000,1);
% 
% magD1 = abs(D1)+0.0002284;
% phaseD1 = angle(D1);
% 
% time1 = 0:length(magD1)/(length(magD1)*6000):length(magD1)/6000 - length(magD1)/(length(magD1)*6000);
% time1 = time1*1.238;
% 
% 
% time = time1+13.3287;%(145388:315021)';%;


% start time
iSSS  = find(abs(time-37.9531)<0.005); iSSS  = iSSS(1); iEEE  = find(abs(time-39)<0.005); iEEE  = iEEE(1);
ibSSS = find(abs(time-42.9531)<0.005); ibSSS = ibSSS(1);ibEEE = find(abs(time-44)<0.005); ibEEE = ibEEE(1);

%x
xSSS = find(abs(time-47.9531)<0.002); xSSS = xSSS(1);xEEE = find(abs(time-51.4225)<0.004); xEEE = xEEE(1);

%y
ySSS = find(abs(time-57.9945)<0.009);   ySSS = ySSS(1);yEEE = find(abs(time-61.0161)<0.009); yEEE = yEEE(1);

%z
zSSS = find(abs(time-68.0237)<0.009); zSSS = zSSS(1);zEEE = find(abs(time-70.034)<0.009); zEEE = zEEE(1);

%xyz back
bSSS1 = find(abs(time-87.9941)<0.009); bSSS1 = bSSS1(1);bEEE1 = find(abs(time-92)<0.009); bEEE1 = bEEE1(1);

%xyz
xyzSSS = find(abs(time-107.99)<0.009);     xyzSSS = xyzSSS(1);xyzEEE = find(abs(time-111.984)<0.009);      xyzEEE = xyzEEE(1);

%xyz back
bSSS = find(abs(time-128)<0.011); bSSS = bSSS(1);   bEEE = find(abs(time-138)<0.009); bEEE = bEEE(1);

timeX = time(xSSS:xEEE); timeY = time(ySSS:yEEE); timeZ = time(zSSS:zEEE); timeXYZ = time(xyzSSS:xyzEEE);

% x then y finally z
[PP1, VV1, AA1] = groundtruth1Dx(timeX-timeX(1)); % move along x
[PP2, VV2, AA2] = groundtruth1Dy(timeY-timeY(1)); % then move along y
[PP3, VV3, AA3] = groundtruth1Dz(timeZ-timeZ(1)); % finally move along z

% xyz together
[PPxyz1, VVxyz1, AAxyz1] = groundtruth1Dx2(timeXYZ-timeXYZ(1)); % move along x y z simulanteously
[PPxyz2, VVxyz2, AAxyz2] = groundtruth1Dy2(timeXYZ-timeXYZ(1)); % move along x y z simulanteously
[PPxyz3, VVxyz3, AAxyz3] = groundtruth1Dz2(timeXYZ-timeXYZ(1)); % move along x y z simulanteously

%
% 3D coordinates
coord3 = zeros(9, length(time));
%
% x
coord3(1,1:iSSS) = 1.03;
coord3(1,iSSS+1:iEEE) = 1.03 + PP3(end)/1.4*(time(iSSS+1:iEEE)-time(iSSS));
coord3(1,iEEE+1:ibSSS) = 1.03 + PP3(end)/1.4;
coord3(1,ibSSS+1:ibEEE) = 1.03 + PP3(end)/1.4 - PP3(end)/1.4*(time(ibSSS+1:ibEEE)-time(ibSSS));

coord3(1,ibEEE:xSSS) = 1.03;
coord3(1,xSSS+1:length(AA1)+xSSS) = PP1 + 1.03;
coord3(1,length(AA1)+xSSS+1:bSSS1) = PP1(end) + 1.03;

coord3(1,bSSS1+1:bSSS1 + length(PPxyz1)) = PP1(end) + 1.03 -  PPxyz1;
coord3(1,bSSS1 + length(PPxyz1)+1:xyzSSS) = 1.03;

coord3(1,xyzSSS+1:xyzSSS+length(PPxyz1)) = PPxyz1 + 1.03;
coord3(1,xyzSSS+1+length(PPxyz1):bSSS) = PPxyz1(end) + 1.03;
coord3(1,bSSS+1:bEEE) = PPxyz1(end) + 1.03 - PP1(end)/10*(time(bSSS+1:bEEE)-time(bSSS));
coord3(1,bEEE+1:end) = 1.03;

%
%y
coord3(2,1:ySSS) = 1.31;

coord3(2,ySSS+1 : length(AA2)+ySSS) = PP2 + 1.31;
coord3(2,length(AA2)+ySSS+1:bSSS1)  = PP2(end) + 1.31;

coord3(2,bSSS1+1:bSSS1 + length(PPxyz2)) = PP2(end) + 1.31 -  PPxyz2;
coord3(2,bSSS1 + length(PPxyz2)+1:xyzSSS) = 1.31;

coord3(2,xyzSSS+1:xyzSSS+length(PPxyz2)) = PPxyz2 + 1.31;
coord3(2,xyzSSS+1+length(PPxyz2):bSSS) = PPxyz2(end) + 1.31;
coord3(2,bSSS+1:bEEE) = PPxyz2(end) + 1.31 - PP2(end)/10*(time(bSSS+1:bEEE)-time(bSSS));
coord3(2,bEEE+1:end) = 1.31;

%z
coord3(3,1:zSSS) = 1.03;
coord3(3,zSSS+1 : zSSS + length(AA3)) = PP3 + 1.03;
coord3(3,length(AA3)+zSSS+1:bSSS1)  = PP3(end) + 1.03;


coord3(3,bSSS1+1:bSSS1 + length(PPxyz3)) = PP3(end) + 1.03 -  PPxyz3;
coord3(3,bSSS1 + length(PPxyz3)+1:xyzSSS) = 1.03;

coord3(3,xyzSSS+1:xyzSSS+length(PPxyz3)) = PPxyz3 + 1.03;
coord3(3,xyzSSS+1+length(PPxyz3):bSSS) = PPxyz3(end) + 1.03;

coord3(3,bSSS+1:bEEE) = PPxyz3(end) + 1.03 - PP3(end)/10*(time(bSSS+1:bEEE)-time(bSSS));
coord3(3,bEEE+1:end) = 1.03;

%===============================================================
%=====================================
% x
coord3(4,iSSS+1:iEEE) = VV3(end)/1.4*(time(iSSS+1:iEEE)-time(iSSS));
coord3(4,iEEE+1:ibSSS) = VV3(end)/1.4;
coord3(4,ibSSS+1:ibEEE) = VV3(end)/1.4 - VV3(end)/1.4*(time(ibSSS+1:ibEEE)-time(ibSSS));

coord3(4,xSSS+1:length(AA1)+xSSS) = VV1 ;
coord3(4,length(AA1)+xSSS+1:bSSS1) = VV1(end);

coord3(4,bSSS1+1:bSSS1 + length(PPxyz1)) = VV1(end) -  VVxyz1;

coord3(4,xyzSSS+1:xyzSSS+length(PPxyz1)) = VVxyz1;
coord3(4,xyzSSS+1+length(PPxyz1):bSSS) = VVxyz1(end);
coord3(4,bSSS+1:bEEE) = VVxyz1(end) - VV1(end)/10*(time(bSSS+1:bEEE)-time(bSSS));
%
%y

coord3(5,ySSS+1 : length(AA2)+ySSS) = VV2;
coord3(5,length(AA2)+ySSS+1:bSSS1)  = VV2(end);

coord3(5,bSSS1+1:bSSS1 + length(PPxyz2)) = VV2(end) -  VVxyz2;

coord3(5,xyzSSS+1:xyzSSS+length(PPxyz2)) = VVxyz2;
coord3(5,xyzSSS+1+length(PPxyz2):bSSS) = VVxyz2(end);
coord3(5,bSSS+1:bEEE) = VVxyz2(end) - VV2(end)/10*(time(bSSS+1:bEEE)-time(bSSS));

%z
coord3(6,zSSS+1 : zSSS + length(AA3)) = VV3;
coord3(6,length(AA3)+zSSS+1:bSSS1)  = VV3(end);


coord3(6,bSSS1+1:bSSS1 + length(PPxyz3)) = VV3(end) -  VVxyz3;

coord3(6,xyzSSS+1:xyzSSS+length(PPxyz3)) = VVxyz3;
coord3(6,xyzSSS+1+length(PPxyz3):bSSS) = VVxyz3(end);
coord3(6,bSSS+1:bEEE) = VVxyz3(end) - VV3(end)/10*(time(bSSS+1:bEEE)-time(bSSS));

%===============================================================
%=====================================
% x
coord3(7,iSSS+1:iEEE) = AA3(end)/1.4*(time(iSSS+1:iEEE)-time(iSSS));
coord3(7,iEEE+1:ibSSS) = AA3(end)/1.4;
coord3(7,ibSSS+1:ibEEE) = AA3(end)/1.4 - AA3(end)/1.4*(time(ibSSS+1:ibEEE)-time(ibSSS));

coord3(7,xSSS+1:length(AA1)+xSSS) = AA1 ;
coord3(7,length(AA1)+xSSS+1:bSSS1) = AA1(end);

coord3(7,bSSS1+1:bSSS1 + length(PPxyz1)) = AA1(end) -  AAxyz1;

coord3(7,xyzSSS+1:xyzSSS+length(PPxyz1)) = AAxyz1;
coord3(7,xyzSSS+1+length(PPxyz1):bSSS) = AAxyz1(end);
coord3(7,bSSS+1:bEEE) = AAxyz1(end) - AA1(end)/10*(time(bSSS+1:bEEE)-time(bSSS));
%
%y

coord3(8,ySSS+1 : length(AA2)+ySSS) = AA2;
coord3(8,length(AA2)+ySSS+1:bSSS1)  = AA2(end);

coord3(8,bSSS1+1:bSSS1 + length(PPxyz2)) = AA2(end) -  AAxyz2;

coord3(8,xyzSSS+1:xyzSSS+length(PPxyz2)) = AAxyz2;
coord3(8,xyzSSS+1+length(PPxyz2):bSSS) = AAxyz2(end);
coord3(8,bSSS+1:bEEE) = AAxyz2(end) - AA2(end)/10*(time(bSSS+1:bEEE)-time(bSSS));

%z
coord3(9,zSSS+1 : zSSS + length(AA3)) = AA3;
coord3(9,length(AA3)+zSSS+1:bSSS1)  = AA3(end);


coord3(9,bSSS1+1:bSSS1 + length(PPxyz3)) = AA3(end) -  AAxyz3;

coord3(9,xyzSSS+1:xyzSSS+length(PPxyz3)) = AAxyz3;
coord3(9,xyzSSS+1+length(PPxyz3):bSSS) = AAxyz3(end);
coord3(9,bSSS+1:bEEE) = AAxyz3(end) - AA3(end)/10*(time(bSSS+1:bEEE)-time(bSSS));
%%
% get z H1 r_sim
% Get RSS and phase from each reader observing the moving tag
z = NaN(3,length(coord3)-1); z_prev = NaN(3,length(coord3)-1);

z_prev(1,:) = coord3(1,1:end-1); z(1,:) = coord3(1,2:end);% x coordinate
z_prev(2,:) = coord3(2,1:end-1); z(2,:) = coord3(2,2:end);% y coordinate
z_prev(3,:) = coord3(3,1:end-1); z(3,:) = coord3(3,2:end);% z coordinate

%
H1 = NaN(1,length(coord3));   H1_ = NaN(1,length(coord3));  r_sim1 = NaN(1,length(coord3)); r_sim1_ = NaN(1,length(coord3)); rdot_sim1 = NaN(1,length(coord3)); rdot_sim1_ = NaN(1,length(coord3)); 
H2 = NaN(1,length(coord3));   H2_ = NaN(1,length(coord3));  r_sim2 = NaN(1,length(coord3)); r_sim2_ = NaN(1,length(coord3)); rdot_sim2 = NaN(1,length(coord3)); rdot_sim2_ = NaN(1,length(coord3)); 
H3 = NaN(1,length(coord3));   H3_ = NaN(1,length(coord3));  r_sim3 = NaN(1,length(coord3)); r_sim3_ = NaN(1,length(coord3)); rdot_sim3 = NaN(1,length(coord3)); rdot_sim3_ = NaN(1,length(coord3)); 
H4 = NaN(1,length(coord3));   H4_ = NaN(1,length(coord3));  r_sim4 = NaN(1,length(coord3)); r_sim4_ = NaN(1,length(coord3)); rdot_sim4 = NaN(1,length(coord3)); rdot_sim4_ = NaN(1,length(coord3));


r_meas1= NaN(1,length(coord3));  rphase1= NaN(1,length(coord3)); phigt1_1= NaN(1,length(coord3)); phigt2_1= NaN(1,length(coord3)); phi1_1 = NaN(1,length(coord3)); phi2_1 = NaN(1,length(coord3));phi3_1 = NaN(1,length(coord3));phi4_1 = NaN(1,length(coord3));
r_meas2= NaN(1,length(coord3));  rphase2= NaN(1,length(coord3)); phigt1_2= NaN(1,length(coord3)); phigt2_2= NaN(1,length(coord3)); phi1_2 = NaN(1,length(coord3)); phi2_2 = NaN(1,length(coord3));phi3_2 = NaN(1,length(coord3));phi4_2 = NaN(1,length(coord3));
r_meas3= NaN(1,length(coord3));  rphase3= NaN(1,length(coord3)); phigt1_3= NaN(1,length(coord3)); phigt2_3= NaN(1,length(coord3)); phi1_3 = NaN(1,length(coord3)); phi2_3 = NaN(1,length(coord3));phi3_3 = NaN(1,length(coord3));phi4_3 = NaN(1,length(coord3));
r_meas4= NaN(1,length(coord3));  rphase4= NaN(1,length(coord3)); phigt1_4= NaN(1,length(coord3)); phigt2_4= NaN(1,length(coord3)); phi1_4 = NaN(1,length(coord3)); phi2_4 = NaN(1,length(coord3));phi3_4 = NaN(1,length(coord3));phi4_4 = NaN(1,length(coord3));

phase1 = NaN(1,length(coord3));
phase2 = NaN(1,length(coord3));
phase3 = NaN(1,length(coord3));
phase4 = NaN(1,length(coord3));

r_sim = NaN(3, length(coord3));  r_meas = NaN(3, length(coord3)); rphase = NaN(3, length(coord3));

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