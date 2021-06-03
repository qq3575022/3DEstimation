clc, clear, close all
[time, coord3] = get3Dcoord();



Gt = 1.462*3/4*2;    % tag's antenna gain
X = 0.35;          % polarization mismatch
M = 4;             % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

x1 = [0,    0,  0.865];  
x2 = [2.29, 0,  1.27];   
x3 = [2.29,2.52, 0.865]; 
x4 = [0, 2.52,  1.27];
%%
radial = NaN(4, length(coord3));
for i = 1:1:length(coord3)
radial(1,i) = sqrt((coord3(1,i) - x1(1))^2+(coord3(2,i) - x1(2))^2+(coord3(3,i) - x1(3))^2);
radial(2,i) = sqrt((coord3(1,i) - x2(1))^2+(coord3(2,i) - x2(2))^2+(coord3(3,i) - x2(3))^2);
radial(3,i) = sqrt((coord3(1,i) - x3(1))^2+(coord3(2,i) - x3(2))^2+(coord3(3,i) - x3(3))^2);
radial(4,i) = sqrt((coord3(1,i) - x4(1))^2+(coord3(2,i) - x4(2))^2+(coord3(3,i) - x4(3))^2);

end
%%
% Parameters of reader
PT = 1;         % reader's transmitted power
GT = 0.01462*0.1;     % reader's trasmitter antenna gain 9.5dBi
GR = 0.01462*0.1;     % reader's receiver   antenna gain 9.5dBi
R = 2;

% Channel noise error covariance
sigma = 0.0000012; 

% phase cconcatenation
% global l1; global l2; global l3; l1 = 0;l2 = 0;l3 = 0;k = 1; 

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

r_sim1 = NaN(1,length(coord3)); r_sim1_ = NaN(1,length(coord3)); rdot_sim1 = NaN(1,length(coord3)); rdot_sim1_ = NaN(1,length(coord3));  phi11 = NaN(1,length(coord3)); phi21 = NaN(1,length(coord3));
r_sim2 = NaN(1,length(coord3)); r_sim2_ = NaN(1,length(coord3)); rdot_sim2 = NaN(1,length(coord3)); rdot_sim2_ = NaN(1,length(coord3));  phi12 = NaN(1,length(coord3)); phi22 = NaN(1,length(coord3));
r_sim3 = NaN(1,length(coord3)); r_sim3_ = NaN(1,length(coord3)); rdot_sim3 = NaN(1,length(coord3)); rdot_sim3_ = NaN(1,length(coord3));  phi13 = NaN(1,length(coord3)); phi23 = NaN(1,length(coord3));
r_sim4 = NaN(1,length(coord3)); r_sim4_ = NaN(1,length(coord3)); rdot_sim4 = NaN(1,length(coord3)); rdot_sim4_ = NaN(1,length(coord3));  phi14 = NaN(1,length(coord3)); phi24 = NaN(1,length(coord3));

%multi = 8*random('Rayleigh',8,[20,1]); gamma = 0.02*rand(20,1) + 0.02*rand(20,1)*i;

multi = randi(10, [20,1]);

for k = 1:1:length(time)-1  
[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1), rdot_sim1(k+1), rdot_sim1_(k+1), phi11(k+1), phi21(k+1)] = noisysim(x1,f1,Gt,M,X,PT,  GT,  0.1*GR,R,sigma,1,k,z,z_prev,time(k+1)-time(k),time(k),multi);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1), rdot_sim2(k+1), rdot_sim2_(k+1), phi12(k+1), phi22(k+1)] = noisysim(x2,f2,Gt,M,X,PT,  GT,  2*GR,R,sigma,2,k,z,z_prev,time(k+1)-time(k),time(k),multi);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1), rdot_sim3(k+1), rdot_sim3_(k+1), phi13(k+1), phi23(k+1)] = noisysim(x3,f3,Gt,M,X,PT,  GT,  6.5*GR,R,sigma,3,k,z,z_prev,time(k+1)-time(k),time(k),multi);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1), rdot_sim4(k+1), rdot_sim4_(k+1), phi14(k+1), phi24(k+1)] = noisysim(x4,f4,Gt,M,X,PT,  GT,  0.04*GR,R,sigma,4,k,z,z_prev,time(k+1)-time(k),time(k),multi);  
end
%%
figure
subplot(411), plot(phi11,'LineWidth',3), hold on, plot(phi21)
subplot(412), plot(phi12,'LineWidth',3), hold on, plot(phi22)
subplot(413), plot(phi13,'LineWidth',3), hold on, plot(phi23)
subplot(414), plot(phi14,'LineWidth',3), hold on, plot(phi24)

%%
phi11_1 = NaN(1,length(coord3));
phi12_1 = NaN(1,length(coord3));
phi13_1 = NaN(1,length(coord3));
phi14_1 = NaN(1,length(coord3));

phi21_2 = NaN(1,length(coord3));
phi22_2 = NaN(1,length(coord3));
phi23_2 = NaN(1,length(coord3));
phi24_2 = NaN(1,length(coord3));

offset1_1 = 0; offset2_1 = 0; offset3_1 = 0; offset4_1 = 0; 
offset1_2 = 0; offset2_2 = 0; offset3_2 = 0; offset4_2 = 0; 

for k = 2:1:length(time)-1  
    [phi11_1(k+1),phi21_2(k+1), offset1_1, offset1_2] = concatenate(phi11, phi21, k, offset1_1, offset1_2);
    [phi12_1(k+1),phi22_2(k+1), offset2_1, offset2_2] = concatenate(phi12, phi22, k, offset2_1, offset2_2);
    [phi13_1(k+1),phi23_2(k+1), offset3_1, offset3_2] = concatenate(phi13, phi23, k, offset3_1, offset3_2);
    [phi14_1(k+1),phi24_2(k+1), offset4_1, offset4_2] = concatenate(phi14, phi24, k, offset4_1, offset4_2);
end

figure
subplot(411), plot(phi11_1), hold on, plot(phi21_2)
subplot(412), plot(phi12_1), hold on, plot(phi22_2)
subplot(413), plot(phi13_1), hold on, plot(phi23_2)
subplot(414), plot(phi14_1), hold on, plot(phi24_2)
%%
for k = 3:1:length(time)-1  
    r_sim1_(k+1) = recover(phi11_1(k), phi21_2(k), f1);
    r_sim2_(k+1) = recover(phi12_1(k), phi22_2(k), f2);
    r_sim3_(k+1) = recover(phi13_1(k), phi23_2(k), f3);
    r_sim4_(k+1) = recover(phi14_1(k), phi24_2(k), f4);
end
%%
figure
subplot(411), plot(r_sim1_);hold on; plot(radial(1,:),'LineWidth',2);
subplot(412), plot(r_sim2_);hold on; plot(radial(2,:),'LineWidth',2);
subplot(413), plot(r_sim3_);hold on; plot(radial(3,:),'LineWidth',2);
subplot(414), plot(r_sim4_);hold on; plot(radial(4,:),'LineWidth',2);
%%
figure
subplot(4,1,1),plot(time, -1.8*r_sim1_+1.67,'LineWidth',3);hold on; plot(time, radial(1,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Recovery from Multi-Path','Ground Truth')
subplot(4,1,2),plot(time, -1.5*r_sim2_+1.73,'LineWidth',3);hold on; plot(time, radial(2,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Recovery from Multi-Path','Ground Truth','location','SouthEast')
subplot(4,1,3),plot(time, -r_sim3_+1.1,'LineWidth',3);hold on; plot(time, radial(3,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Recovery from Multi-Path','Ground Truth','location','SouthEast')
subplot(4,1,4),plot(time, -1.4*r_sim4_+0.8,'LineWidth',3);hold on; plot(time, radial(4,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Recovery from Multi-Path','Ground Truth')

%%
% figure
% subplot(4,1,1),plot(time, H1,'LineWidth',2);hold on; plot(time, H1_,'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,2),plot(time, H2,'LineWidth',2);hold on; plot(time, H2_,'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,3),plot(time, H3,'LineWidth',2);hold on; plot(time, H3_,'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,4),plot(time, H4,'LineWidth',2);hold on; plot(time, H4_,'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% 
% %%
% figure
% subplot(4,1,1),plot(time, r_sim1,'LineWidth',2);hold on; plot(time, r_sim1_,'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,2),plot(time, r_sim2,'LineWidth',2);hold on; plot(time, r_sim2_,'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,3),plot(time, r_sim3,'LineWidth',2);hold on; plot(time, r_sim3_,'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,4),plot(time, r_sim4,'LineWidth',2);hold on; plot(time, r_sim4_,'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
%%
figure
subplot(4,1,1),plot(time, 0.26*r_sim1_-1.92,'LineWidth',3);hold on; plot(time, radial(1,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Recovery from Multi-Path','Ground Truth')
subplot(4,1,2),plot(time, 0.18*r_sim2_-1.5,'LineWidth',3);hold on; plot(time, radial(2,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Recovery from Multi-Path','Ground Truth','location','SouthEast')
subplot(4,1,3),plot(time, -0.53*r_sim3_+11.5,'LineWidth',3);hold on; plot(time, radial(3,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Recovery from Multi-Path','Ground Truth','location','SouthEast')
subplot(4,1,4),plot(time, 0.21*r_sim4_-1.25,'LineWidth',3);hold on; plot(time, radial(4,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Recovery from Multi-Path','Ground Truth')

%%
figure
subplot(4,1,1),plot(time, rdot_sim1,'LineWidth',3);hold on; plot(time, radial(1,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,2),plot(time, rdot_sim2,'LineWidth',3);hold on; plot(time, radial(2,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,3),plot(time, rdot_sim3,'LineWidth',3);hold on; plot(time, radial(3,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,4),plot(time, rdot_sim4,'LineWidth',3);hold on; plot(time, radial(4,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')

%%
% Reader 1
D0_1 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader1_6.bin', 1000000000, 1);
D1 = D0_1;%(85000:445000,1);

magD1 = abs(D1);
phaseD1 = angle(D1);

time1 = 0:length(magD1)/(length(magD1)*6000):length(magD1)/6000 - length(magD1)/(length(magD1)*6000);
time1 = time1*1.238;

[time2_1, magD3_1, magD4_1, magD5_1, magD6_1, magD7_1, phaseD2_1] = filterAvg(magD1, phaseD1);


% Reader 2
D02 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader2_6.bin', 1000000000, 1);
D2 = D02;%(152000:512000,1);

magD2 = abs(D2);
phaseD2 = angle(D2);

time2 = 0:length(magD2)/(length(magD2)*6000):length(magD2)/6000 - length(magD2)/(length(magD2)*6000);
time2 = time2*1.238;

[time2_2, magD3_2, magD4_2, magD5_2, magD6_2, magD7_2, phaseD2_2] = filterAvg(magD2, phaseD2);

% Reader 3
D03 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader3_6.bin', 1000000000, 1);
D3 = D03(31558:end,1);

magD3 = abs(D3);
phaseD3 = angle(D3);


time3 = 0:length(magD3)/(length(magD3)*12000):length(magD3)/12000 - length(magD3)/(length(magD3)*12000);
time3 = time3*1.217;

[time2_3, magD3_3, magD4_3, magD5_3, magD6_3, magD7_3, phaseD2_3] = filterAvg2(magD3, phaseD3);

% Reader 4
D04 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader4_6.bin', 1000000000, 1);
D4 = D04;%(150000:510000,1);

magD4 = abs(D4);
phaseD4 = angle(D4);

time4 = 0:length(magD4)/(length(magD4)*6000):length(magD4)/6000 - length(magD4)/(length(magD4)*6000);
time4 = time4*1.238;

%save('RF.mat','t1','mag1','phase1','t2','mag2','phase2','t3','mag3','phase3','t4','mag4','phase4');

%[time2_4, magD3_4, magD4_4, magD5_4, magD6_4, magD7_4, phaseD2_4] = filterAvg(magD4, phaseD4);

%
mag1   = magD1(145388:315021);     mag2 = magD2(145388:315021);     mag3 = magD3(295791:640912);     mag4 = magD4(145388:315021);
phase1 = phaseD1(145388:315021); phase2 = phaseD2(145388:315021); phase3 = phaseD3(295791:640912); phase4 = phaseD4(145388:315021);
t1 = time1(145388:315021);           t2 = time2(145388:315021);       t3 = time3(295791:640912);      t4 = time4(145388:315021); 

% 
% figure
% subplot(411), plot(time1, phaseD1);
% subplot(412), plot(time2, phaseD2);
% subplot(413), plot(time3, phaseD3);
% subplot(414), plot(time4, phaseD4);
% 
% figure
% subplot(411), plot(time1, magD1)
% subplot(412), plot(time2, magD2)
% subplot(413), plot(time3, magD3)
% subplot(414), plot(time4, magD4)
%
magD12 = NaN(1,length(time)); magD22 = NaN(1,length(time)); magD32 = NaN(1,length(time)); magD42 = NaN(1,length(time));
magD12(1) = magD1(1); magD22(1) = magD1(1);  magD32(1) = magD1(1);  magD42(1) = magD1(1); 
indexT = 2;

for indexMag = 1:1:length(time1)
%     timeT = time(indexT)
%     time1T = time1(indexMag)+13.3287
    if abs(time1(indexMag)+13.3287 - time(indexT)) < 0.0002 || ((time1(indexMag)+13.3287) > time(indexT))
        %indexT
        if indexT  <= length(time)
            magD12(indexT) = magD1(indexMag);
            magD22(indexT) = magD2(indexMag);
            
            if indexMag < length(time4)
                magD42(indexT) = magD4(indexMag);
            end
        end
        
        indexT = indexT + 1;
    end
end

indexT = 1;

for indexMag = 1:1:length(time3)
    if abs(time3(indexMag)+13.3287 - time(indexT)) < 0.0002 || ((time3(indexMag)+13.3287) > time(indexT))
        
        if indexT  <= length(time)
            magD32(indexT) = magD3(indexMag);
        end
        
        indexT = indexT + 1;
    end
end
% Measurement Magnitude
% figure
% subplot(411), plot(time, magD12);
% subplot(412), plot(time, magD22);
% subplot(413), plot(time, magD32);
% subplot(414), plot(time, magD42);
%
%

% figure 
% subplot(411), plot(time, magD12, 'LineWidth', 2), title('Reader 1 Magnitude'),xlim([35, 155]);
% hold on
% plot(time, H1, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
% hold on
% legend('Measurement','Simulation','location','SouthEast')
% grid on 
% grid minor
% 
% subplot(412), plot(time, magD22, 'LineWidth', 2), title('Reader 2 Magnitude'),xlim([35, 155]);
% hold on
% plot(time, H2, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
% hold on
% legend('Measurement','Simulation','location','NorthEast')
% grid on 
% grid minor
% 
% subplot(413), plot(time, magD32, 'LineWidth', 2), title('Reader 3 Magnitude'),xlim([35, 155]);
% hold on
% plot(time, H3, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
% hold on
% legend('Measurement','Simulation','location','NorthEast')
% grid on 
% grid minor
% 
% subplot(414), plot(time, magD42, 'LineWidth', 2), title('Reader 4 Magnitude'),xlim([35, 155]);
% hold on
% plot(time, H4, 'LineWidth', 2),xlim([35, 155]);
% hold on
% legend('Measurement','Simulation','location','NorthEast'), ylabel('Magnitude [V]');
% grid on 
% grid minor

%

[-mean(H1(2:124840) - magD12(2:124840)), rms(H1(2:124840) - magD12(2:124840)), sqrt(rms(H1(2:124840) - magD12(2:124840)))]
[-mean(H2(2:124840) - magD22(2:124840)), rms(H2(2:124840) - magD22(2:124840)), sqrt(rms(H2(2:124840) - magD22(2:124840)))]
[-mean(H3(2:124840) - magD32(2:124840)), rms(H3(2:124840) - magD32(2:124840)), sqrt(rms(H3(2:124840) - magD32(2:124840)))]
[-mean(H4(2:124840) - magD42(2:124840)), rms(H4(2:124840) - magD42(2:124840)), sqrt(rms(H4(2:124840) - magD42(2:124840)))]


[-mean(H1_(2:124840) - magD12(2:124840)), rms(H1_(2:124840) - magD12(2:124840)), sqrt(rms(H1_(2:124840) - magD12(2:124840)))]
[-mean(H2_(2:124840) - magD22(2:124840)), rms(H2_(2:124840) - magD22(2:124840)), sqrt(rms(H2_(2:124840) - magD22(2:124840)))]
[-mean(H3_(2:124840) - magD32(2:124840)), rms(H3_(2:124840) - magD32(2:124840)), sqrt(rms(H3_(2:124840) - magD32(2:124840)))]
[-mean(H4_(2:124840) - magD42(2:124840)), rms(H4_(2:124840) - magD42(2:124840)), sqrt(rms(H4_(2:124840) - magD42(2:124840)))]

[-mean(H1_(2:124840) - H1(2:124840)), rms(H1_(2:124840) - H1(2:124840)), sqrt(rms(H1_(2:124840) - H1(2:124840)))]
[-mean(H2_(2:124840) - H2(2:124840)), rms(H2_(2:124840) - H2(2:124840)), sqrt(rms(H2_(2:124840) - H2(2:124840)))]
[-mean(H3_(2:124840) - H3(2:124840)), rms(H3_(2:124840) - H3(2:124840)), sqrt(rms(H3_(2:124840) - H3(2:124840)))]
[-mean(H4_(2:124840) - H4(2:124840)), rms(H4_(2:124840) - H4(2:124840)), sqrt(rms(H4_(2:124840) - H4(2:124840)))]


%% Get K Value
% m41 = getm(H1);
% m42 = getm(H2);
% m43 = getm(H3);
% m44 = getm(H4);
% 
% figure
% subplot(411), plot(time(1:end-10000), medfilt1(m41, 10000),'LineWidth', 2);title('K Value for Reader 1 in Simulation');ylabel('K Value'); grid on; grid minor;ylim([0,100]);xlim([5,140])
% subplot(412), plot(time(1:end-10000), medfilt1(m42, 10000),'LineWidth', 2);title('K Value for Reader 2 in Simulation');ylabel('K Value'); grid on; grid minor;ylim([70,170]);xlim([5,140])
% subplot(413), plot(time(1:end-10000), medfilt1(m43, 10000),'LineWidth', 2);title('K Value for Reader 3 in Simulation');ylabel('K Value'); grid on; grid minor;ylim([0,600]);xlim([5,140])
% subplot(414), plot(time(1:end-10000), medfilt1(m44, 10000),'LineWidth', 2);title('K Value for Reader 4 in Simulation');ylabel('K Value'); grid on; grid minor;ylim([0,100]);xlim([5,140]);xlabel('time [s]');

%%
% figure
% subplot(411), plot(time1+13.3287, phaseD1,'LineWidth',2), xlim([14, 145]), hold on, plot(time, rdot_sim1+1.57,'LineWidth',1), xlim([14, 145]), legend('Measurement','Simulation','location','BestOutside'), title('Phase from Reader $\#1$ in 3D','interpreter','latex'), ylabel('Phase [rad]'), xlabel('time [s]')
% subplot(412), plot(time2+13.3287, phaseD2,'LineWidth',2), xlim([14, 145]), hold on, plot(time, rdot_sim2+1.57,'LineWidth',1), xlim([14, 145]), legend('Measurement','Simulation','location','BestOutside'), title('Phase from Reader $\#2$ in 3D','interpreter','latex'), ylabel('Phase [rad]'), xlabel('time [s]')
% subplot(413), plot(time3+13.6187, phaseD3,'LineWidth',2), xlim([14, 145]), hold on, plot(time, rdot_sim3+1.56,'LineWidth',1), xlim([14, 145]), legend('Measurement','Simulation','location','BestOutside'), title('Phase from Reader $\#3$ in 3D','interpreter','latex'), ylabel('Phase [rad]'), xlabel('time [s]')
% subplot(414), plot(time4+13.3287, phaseD4,'LineWidth',2), xlim([14, 145]), hold on, plot(time, rdot_sim4+1.55,'LineWidth',1), xlim([14, 145]), legend('Measurement','Simulation','location','BestOutside'), title('Phase from Reader $\#4$ in 3D','interpreter','latex'), ylabel('Phase [rad]'), xlabel('time [s]')

%
% m21 = NaN(1, length(D1) - 10000);
% m31 = NaN(1, length(D2) - 10000);
% m41 = NaN(1, length(D3) - 10000);
% m51 = NaN(1, length(D4) - 10000);
% 
% %m4 = NaN(1, length(D) - 10000);
% 
% for i = 1:1:length(D1) - 10000
%     moment11 = abs(mean(D1(i:i+10000)));
%     moment12 = abs(mean(D1(i:i+10000).^2));
%     
%     m21(i) = (-6*moment11.^2+4*moment12-sqrt(12*moment11.^4-8*moment11.^2.*moment12))./(moment11.^2-moment12);%abs(moment1^2/(moment2 - moment1^2));
% end
% 

% for i = 1:1:length(D2) - 10000
%     moment21 = abs(mean(D2(i:i+10000)));
%     moment22 = abs(mean(D2(i:i+10000).^2));
%     
%     m31(i) = (moment21 - moment22.^2)/4;%abs(moment1^2/(moment2 - moment1^2))
% end
% 
% for i = 1:1:length(D3) - 10000
%     moment31 = abs(mean(D3(i:i+10000)));
%     moment32 = abs(mean(D3(i:i+10000).^2));
% 
%     m41(i) = (moment31 - moment32.^2)/4;%abs(moment1^2/(moment2 - moment1^2));
%         %m4(i) = abs((-2*moment4^2 + moment8 - moment4*sqrt(abs(2*moment4^2 - moment8)))/(moment4^2 - moment8));
% end
% 
% for i = 1:1:length(D4) - 10000
%     moment41 = abs(mean(D4(i:i+10000)));
%     moment42 = abs(mean(D4(i:i+10000).^2));
% 
%     m51(i) = (moment41 - moment42.^2)/4;%abs(moment1^2/(moment2 - moment1^2));
% end

%
% figure
% subplot(411), plot(time1(1:end-10000), medfilt1(m21, 10000),'LineWidth', 2); title('K Value for Reader 1');ylabel('K Value');%legend('m2', 'm4');
% subplot(412), plot(time2(1:end-10000), medfilt1(m31, 10000),'LineWidth', 2); title('K Value for Reader 2');ylabel('K Value');%legend('m2', 'm4');
% subplot(413), plot(time3(1:end-10000), medfilt1(m41, 10000),'LineWidth', 2); title('K Value for Reader 3');ylabel('K Value');%legend('m2', 'm4');
% subplot(414), plot(time4(1:end-10000), medfilt1(m51, 10000),'LineWidth', 2); title('K Value for Reader 4');ylabel('K Value');xlabel('time [s]')%legend('m2', 'm4');

