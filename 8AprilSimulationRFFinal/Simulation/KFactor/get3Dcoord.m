function [time, coord3, z, z_prev, H1, H2, H3, H4, H1_, H2_, H3_, H4_, r_sim1, r_sim2, r_sim3, r_sim4, r_sim1_, r_sim2_, r_sim3_, r_sim4_, rdot_sim1, rdot_sim2, rdot_sim3, rdot_sim4, rdot_sim1_, rdot_sim2_, rdot_sim3_, rdot_sim4_, diff1, diff2, diff3, diff4] = get3Dcoord()

data=readtable('2.csv','Delimiter', ',');  g=9.7953;

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)
sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 
sensor_gyro_check = ismember(data.Var2,'GYRO_UN');  find_gyro_data    = find(sensor_gyro_check == 1);  
sensor_mag_check  = ismember(data.Var2,'MAG_UN');   find_mag_data     = find(sensor_mag_check == 1);

time   = table2array(data(:,1));  data_x = table2array(data(:,3));  data_y = table2array(data(:,4));  data_z = table2array(data(:,5));

acc_time   = time(find_acc_data);  acc_time   = (acc_time - acc_time(1))/1000000000;                
gyro_time  = time(find_gyro_data); gyro_time  = (gyro_time - gyro_time(1))/1000000000;
mag_time   = time(find_mag_data);  mag_time   = (mag_time - mag_time(1))/1000000000;

time  = unique(sort([acc_time; gyro_time; mag_time]),'rows');

% start time
iSSS = find(abs(time-37.9531)<0.002); iSSS = iSSS(1);iEEE = find(abs(time-39)<0.004); iEEE = iEEE(1);
ibSSS = find(abs(time-42.9531)<0.002); ibSSS = ibSSS(1);ibEEE = find(abs(time-44)<0.004); ibEEE = ibEEE(1);


%x
xS = find(abs(acc_time-47.9531)<0.002); xS = xS(1);xE = find(abs(acc_time-51.4225)<0.001); xE = xE(1);
xxS = find(abs(gyro_time-47.9531)<0.002); xxS = xxS(1);xxE = find(abs(gyro_time-51.4225)<0.001); xxE = xxE(1);
xSS = find(abs(mag_time-47.9531)<0.005); xSS = xSS(1);xEE = find(abs(mag_time-51.4225)<0.012); xEE = xEE(1);
xSSS = find(abs(time-47.9531)<0.002); xSSS = xSSS(1);xEEE = find(abs(time-51.4225)<0.004); xEEE = xEEE(1);

%y
yS = find(abs(acc_time-57.9945)<0.002); yS = yS(1);yE = find(abs(acc_time-61.0161)<0.002); yE = yE(1);
yyS = find(abs(gyro_time-57.9945)<0.002);  yyS = yyS(1);yyE = find(abs(gyro_time-61.0161)<0.002); yyE = yyE(1);
ySS = find(abs(mag_time-57.9945)<0.009);   ySS = ySS(1);yEE = find(abs(mag_time-61.0161)<0.012); yEE = yEE(1);
ySSS = find(abs(time-57.9945)<0.009);   ySSS = ySSS(1);yEEE = find(abs(time-61.0161)<0.009); yEEE = yEEE(1);

%z
zS = find(abs(acc_time-68.0237)<0.002); zS = zS(1);zE = find(abs(acc_time-70.034)<0.002); zE = zE(1);
zzS = find(abs(gyro_time-68.0237)<0.002); zzS = zzS(1);zzE = find(abs(gyro_time-70.034)<0.002); zzE = zzE(1);
zSS = find(abs(mag_time-68.0237)<0.009); zSS = zSS(1);zEE = find(abs(mag_time-70.034)<0.009); zEE = zEE(1);
zSSS = find(abs(time-68.0237)<0.009); zSSS = zSSS(1);zEEE = find(abs(time-70.034)<0.009); zEEE = zEEE(1);

%xyz back
bS1 = find(abs(acc_time-87.9941)<0.002); bS1 = bS1(1);bE1 = find(abs(acc_time-92)<0.002); bE1 = bE1(1);
bbS1 = find(abs(gyro_time-87.9941)<0.002); bbS1 = bbS1(1);bbE1 = find(abs(gyro_time-92)<0.002); bbE1 = bbE1(1);
bSS1 = find(abs(mag_time-87.9941)<0.012); bSS1 = bSS1(1);bEE1 = find(abs(mag_time-92)<0.012); bEE1 = bEE1(1);
bSSS1 = find(abs(time-87.9941)<0.009); bSSS1 = bSSS1(1);bEEE1 = find(abs(time-92)<0.009); bEEE1 = bEEE1(1);

%xyz
xyzS = find(abs(acc_time-107.99)<0.002);   xyzS = xyzS(1);    xyzE = find(abs(acc_time-111.984)<0.002);    xyzE = xyzE(1);
xyzzS = find(abs(gyro_time-107.99)<0.002); xyzzS = xyzzS(1);  xyzzE = find(abs(gyro_time-111.984)<0.002);  xyzzE = xyzzE(1);
xyzSS = find(abs(mag_time-107.99)<0.012);  xyzSS = xyzSS(1);  xyzEE = find(abs(mag_time-111.984)<0.009);   xyzEE = xyzEE(1);
xyzSSS = find(abs(time-107.99)<0.009);     xyzSSS = xyzSSS(1);xyzEEE = find(abs(time-111.984)<0.009);      xyzEEE = xyzEEE(1);

%xyz back
bS = find(abs(acc_time-128)<0.002); bS = bS(1);     bE = find(abs(acc_time-138)<0.002); bE = bE(1);
bbS = find(abs(gyro_time-128)<0.002); bbS = bbS(1); bbE = find(abs(gyro_time-138)<0.002); bbE = bbE(1);
bSS = find(abs(mag_time-128)<0.012); bSS = bSS(1);  bEE = find(abs(mag_time-138)<0.012); bEE = bEE(1);
bSSS = find(abs(time-128)<0.009); bSSS = bSSS(1);   bEEE = find(abs(time-138)<0.009); bEEE = bEEE(1);


% ==========   Raw Acceleration Data  ========== 
% ----acc_data_x----
% ----acc_data_y----
% ----acc_data_z----
acc_data_x  = data_x(find_acc_data);   acc_data_y = data_y(find_acc_data)-g; acc_data_z = data_z(find_acc_data);
gyro_data_x = data_x(find_gyro_data);  gyro_data_y = data_y(find_gyro_data); gyro_data_z = data_z(find_gyro_data);
mag_data_x  = data_x(find_mag_data);   mag_data_y = data_y(find_mag_data);   mag_data_z = data_z(find_mag_data);

% 
acc_data_x = medfilt1(acc_data_x,80);
acc_data_y = medfilt1(acc_data_y,80);
acc_data_z = medfilt1(acc_data_z,80);
% 
%x
accxX = acc_data_x(xS:xE) - mean(acc_data_x(xS:xE)); accyX = acc_data_y(xS:xE) - mean(acc_data_y(xS:xE)); acczX = acc_data_z(xS:xE) - mean(acc_data_z(xS:xE));
accTX = acc_time(xS:xE);

gyroxX = gyro_data_x(xxS:xxE) - mean(gyro_data_x(xxS:xxE)); gyroyX = gyro_data_y(xxS:xxE) - mean(gyro_data_y(xxS:xxE)); gyrozX = gyro_data_z(xxS:xxE) - mean(gyro_data_z(xxS:xxE));
gyroTX = gyro_time(xxS:xxE);

magxX = mag_data_x(xSS:xEE) - mean(mag_data_x(xSS:xEE)); magyX = mag_data_y(xSS:xEE) - mean(mag_data_y(xSS:xEE)); magzX = mag_data_z(xSS:xEE) - mean(mag_data_z(xSS:xEE));
magTX = mag_time(xSS:xEE);
% timeX
timeX = unique(sort([accTX; gyroTX; magTX]),'rows');

%y
accxY = acc_data_x(yS:yE) - mean(acc_data_x(yS:yE)); accyY = acc_data_y(yS:yE) - mean(acc_data_y(yS:yE)); acczY = acc_data_z(yS:yE) - mean(acc_data_z(yS:yE));
accTY = acc_time(yS:yE);

gyroxY = gyro_data_x(yyS:yyE) - mean(gyro_data_x(yyS:yyE)); gyroyY = gyro_data_y(yyS:yyE) - mean(gyro_data_y(yyS:yyE)); gyrozY = gyro_data_z(yyS:yyE) - mean(gyro_data_z(yyS:yyE));
gyroTY = gyro_time(yyS:yyE);

magxY = mag_data_x(ySS:yEE) - mean(mag_data_x(ySS:yEE)); magyY = mag_data_y(ySS:yEE) - mean(mag_data_y(ySS:yEE)); magzY = mag_data_z(ySS:yEE) - mean(mag_data_z(ySS:yEE));
magTY = mag_time(ySS:yEE);
% timeY
timeY = unique(sort([accTY; gyroTY; magTY]),'rows');

%z
accxZ = acc_data_x(zS:zE) - mean(acc_data_x(zS:zE)); accyZ = acc_data_y(zS:zE) - mean(acc_data_y(zS:zE)); acczZ = acc_data_z(zS:zE) - mean(acc_data_z(zS:zE));
accTZ = acc_time(zS:zE);

gyroxZ = gyro_data_x(zzS:zzE) - mean(gyro_data_x(zzS:zzE)); gyroyZ = gyro_data_y(zzS:zzE) - mean(gyro_data_y(zzS:zzE)); gyrozZ = gyro_data_z(zzS:zzE) - mean(gyro_data_z(zzS:zzE));
gyroTZ = gyro_time(zzS:zzE);

magxZ = mag_data_x(zSS:zEE) - mean(mag_data_x(zSS:zEE)); magyZ = mag_data_y(zSS:zEE) - mean(mag_data_y(zSS:zEE)); magzZ = mag_data_z(zSS:zEE) - mean(mag_data_z(zSS:zEE));
magTZ = mag_time(zSS:zEE);
% timeZ
timeZ = unique(sort([accTZ; gyroTZ; magTZ]),'rows');

%xyz back
accxbZ = acc_data_x(bS1:bE1) - mean(acc_data_x(bS1:bE1)); accybZ = acc_data_y(bS1:bE1) - mean(acc_data_y(bS1:bE1)); acczbZ = acc_data_z(bS1:bE1) - mean(acc_data_z(bS1:bE1));
accTbZ = acc_time(bS1:bE1);

gyroxbZ = gyro_data_x(bbS1:bbE1) - mean(gyro_data_x(bbS1:bbE1)); gyroybZ = gyro_data_y(bbS1:bbE1) - mean(gyro_data_y(bbS1:bbE1)); gyrozbZ = gyro_data_z(bbS1:bbE1) - mean(gyro_data_z(bbS1:bbE1));
gyroTbZ = gyro_time(bbS1:bbE1);

magxbZ = mag_data_x(bSS1:bEE1) - mean(mag_data_x(bSS1:bEE1)); magybZ = mag_data_y(bSS1:bEE1) - mean(mag_data_y(bSS1:bEE1)); magzbZ = mag_data_z(bSS1:bEE1) - mean(mag_data_z(bSS1:bEE1));
magTbZ = mag_time(bSS1:bEE1);

% timeZ
timebZ = unique(sort([accTbZ; gyroTbZ; magTbZ]),'rows');

% xyz
accxZ = acc_data_x(xyzS:xyzE) - mean(acc_data_x(xyzS:xyzE)); accyZ = acc_data_y(xyzS:xyzE) - mean(acc_data_y(xyzS:xyzE)); acczZ = acc_data_z(xyzS:xyzE) - mean(acc_data_z(xyzS:xyzE));
accTXYZ = acc_time(xyzS:xyzE);

gyroxZ = gyro_data_x(xyzzS:xyzzE) - mean(gyro_data_x(xyzzS:xyzzE)); gyroyZ = gyro_data_y(xyzzS:xyzzE) - mean(gyro_data_y(xyzzS:xyzzE)); gyrozZ = gyro_data_z(xyzzS:xyzzE) - mean(gyro_data_z(xyzzS:xyzzE));
gyroTXYZ = gyro_time(xyzzS:xyzzE);

magxZ = mag_data_x(xyzSS:xyzEE) - mean(mag_data_x(xyzSS:xyzEE)); magyZ = mag_data_y(xyzSS:xyzEE) - mean(mag_data_y(xyzSS:xyzEE)); magzZ = mag_data_z(xyzSS:xyzEE) - mean(mag_data_z(xyzSS:xyzEE));
magTXYZ = mag_time(xyzSS:xyzEE);

% timeXYZ
timeXYZ = unique(sort([accTXYZ; gyroTXYZ; magTXYZ]),'rows');

% Acceleration

% x then y finally z
[P1, V1, A1] = groundtruth1Dx(accTX - accTX(1));
[P2, V2, A2] = groundtruth1Dy(accTY - accTY(1));
[P3, V3, A3] = groundtruth1Dz(accTZ - accTZ(1));

% xyz together
[Pxyz1, Vxyz1, Axyz1] = groundtruth1Dx2(accTXYZ-accTXYZ(1));
[Pxyz2, Vxyz2, Axyz2] = groundtruth1Dy2(accTXYZ-accTXYZ(1));
[Pxyz3, Vxyz3, Axyz3] = groundtruth1Dz2(accTXYZ-accTXYZ(1));

%
AA =  zeros(3, length(acc_time));

% x
AA(1,xS+1:length(A1) + xS) = A1; AA(1,xyzS+1:length(Axyz1) + xyzS) = Axyz1; %AA(1,bS1+1:length(Axyz1) + bS1) = -Axyz1;%A3(1,xyzSback+1:length(AA1) + xyzSback) = -AA1; A3(1,xyzS+1:length(AA1) + xyzS) = AA1;
%y
AA(2,yS+1:length(A2) + yS) = A2; AA(2,xyzS+1:length(Axyz2) + xyzS) = Axyz2; %AA(2,bS1+1:length(Axyz2) + bS1) = -Axyz2;%A3(2,xyzSback+1:length(AA2) + xyzSback) = -AA2; A3(2,xyzS+1:length(AA2) + xyzS) = AA2;
%z
AA(3,zS+1:length(A3) + zS) = A3; AA(3,xyzS+1:length(Axyz3) + xyzS) = Axyz3; %AA(3,bS1+1:length(Axyz3) + bS1) = -Axyz3;%A3(3,xyzSback+1:length(AA3) + xyzSback) = -AA3; A3(3,xyzS+1:length(AA3) + xyzS) = AA3;
%time = unique(sort([accT; gyroT; magT]),'rows');

% All Measurement

% x then y finally z
[PP1, VV1, AA1] = groundtruth1Dx(timeX-timeX(1));
[PP2, VV2, AA2] = groundtruth1Dy(timeY-timeY(1));
[PP3, VV3, AA3] = groundtruth1Dz(timeZ-timeZ(1));

% xyz together
[PPxyz1, VVxyz1, AAxyz1] = groundtruth1Dx2(timeXYZ-timeXYZ(1));
[PPxyz2, VVxyz2, AAxyz2] = groundtruth1Dy2(timeXYZ-timeXYZ(1));
[PPxyz3, VVxyz3, AAxyz3] = groundtruth1Dz2(timeXYZ-timeXYZ(1));

%
% 3D coordinates
coord3 = zeros(3, length(time));
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


%% get z H1 r_sim

% Get RSS and phase from each reader observing the moving tag
z = NaN(3,length(coord3)-1); z_prev = NaN(3,length(coord3)-1);

z_prev(1,:) = coord3(1,1:end-1); z(1,:) = coord3(1,2:end);% x coordinate
z_prev(2,:) = coord3(2,1:end-1); z(2,:) = coord3(2,2:end);% y coordinate
z_prev(3,:) = coord3(3,1:end-1); z(3,:) = coord3(3,2:end);% z coordinate

%
H1 = NaN(1,length(coord3));   H1_ = NaN(1,length(coord3));  %phi1 = NaN(1,length(coord3));       phi_mu1 = NaN(1,length(coord3));
H2 = NaN(1,length(coord3));   H2_ = NaN(1,length(coord3));   %phi2 = NaN(1,length(coord3));       phi_mu2 = NaN(1,length(coord3));
H3 = NaN(1,length(coord3));   H3_ = NaN(1,length(coord3)); %phi3 = NaN(1,length(coord3));       phi_mu3 = NaN(1,length(coord3));
H4 = NaN(1,length(coord3));   H4_ = NaN(1,length(coord3));  %phi4 = NaN(1,length(coord3));       phi_mu4 = NaN(1,length(coord3));

r_sim1 = NaN(1,length(coord3)); r_sim1_ = NaN(1,length(coord3)); rdot_sim1 = NaN(1,length(coord3)); rdot_sim1_ = NaN(1,length(coord3));  diff1 = NaN(1,length(coord3));
r_sim2 = NaN(1,length(coord3)); r_sim2_ = NaN(1,length(coord3)); rdot_sim2 = NaN(1,length(coord3)); rdot_sim2_ = NaN(1,length(coord3));  diff2 = NaN(1,length(coord3));
r_sim3 = NaN(1,length(coord3)); r_sim3_ = NaN(1,length(coord3)); rdot_sim3 = NaN(1,length(coord3)); rdot_sim3_ = NaN(1,length(coord3));  diff3 = NaN(1,length(coord3));
r_sim4 = NaN(1,length(coord3)); r_sim4_ = NaN(1,length(coord3)); rdot_sim4 = NaN(1,length(coord3)); rdot_sim4_ = NaN(1,length(coord3));  diff4 = NaN(1,length(coord3));


end