function [coord3, radial, z, z_prev, H1, H2, H3, H4, H1_, H2_, H3_, H4_, r_sim, r_sim1, r_sim2, r_sim3, r_sim4, r_sim1_, r_sim2_, r_sim3_, r_sim4_, r_meas, r_meas1, r_meas2, r_meas3, r_meas4, rphase, rphase1, rphase2, rphase3, rphase4, rdot_sim1, rdot_sim2, rdot_sim3, rdot_sim4, rdot_sim1_, rdot_sim2_, rdot_sim3_, rdot_sim4_, phase1, phase2, phase3, phase4, phigt1_1, phigt1_2, phigt1_3, phigt1_4, phigt2_1, phigt2_2, phigt2_3, phigt2_4, phi1_1, phi1_2, phi1_3, phi1_4, phi2_1, phi2_2, phi2_3, phi2_4, phi3_1, phi3_2, phi3_3, phi3_4, phi4_1, phi4_2, phi4_3, phi4_4] = get3Dcoord(x1, x2, x3, x4, time)

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