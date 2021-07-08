clc, clear, close all
[time, coord3, z, z_prev, H1, H2, H3, H4, H1_, H2_, H3_, H4_, r_sim1, r_sim2, r_sim3, r_sim4, r_sim1_, r_sim2_, r_sim3_, r_sim4_, rdot_sim1, rdot_sim2, rdot_sim3, rdot_sim4, rdot_sim1_, rdot_sim2_, rdot_sim3_, rdot_sim4_, phase1, phase2, phase3, phase4, phi1, phi2, phi3, phi4 , phi21, phi22, phi23, phi24, phi31, phi32, phi33, phi34 , phi41, phi42, phi43, phi44] = get3Dcoord();

x1 = [0,    0,  0.865];  
x2 = [2.29, 0,  1.27];   
x3 = [2.29,2.52, 0.865]; 
x4 = [0, 2.52,  1.27];

radial = NaN(4, length(coord3));

for i = 1:1:length(coord3)
radial(1,i) = sqrt((coord3(1,i) - x1(1))^2+(coord3(2,i) - x1(2))^2+(coord3(3,i) - x1(3))^2);
radial(2,i) = sqrt((coord3(1,i) - x2(1))^2+(coord3(2,i) - x2(2))^2+(coord3(3,i) - x2(3))^2);
radial(3,i) = sqrt((coord3(1,i) - x3(1))^2+(coord3(2,i) - x3(2))^2+(coord3(3,i) - x3(3))^2);
radial(4,i) = sqrt((coord3(1,i) - x4(1))^2+(coord3(2,i) - x4(2))^2+(coord3(3,i) - x4(3))^2);

end
%
% Parameters of reader
Gt = 1/75*sqrt(1.462*3/4);    % tag's antenna gain
X = 0.85;                     % polarization mismatch
M = sqrt(4);                  % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

% Parameters of reader
PT = 1;                      % reader's transmitted power
GT = sqrt(1.462);            % reader's trasmitter antenna gain 9.5dBi
GR = sqrt(1.462);            % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.0000012; 

% phase cconcatenation
% global l1; global l2; global l3; l1 = 0;l2 = 0;l3 = 0;k = 1; 


offset11 = 0; offset12 = 0; offset13 = 0; offset14 = 0;  offset21 = 0; offset22 = 0; offset23 = 0; offset24 = 0;
offset31 = 0; offset32 = 0; offset33 = 0; offset34 = 0;  offset41 = 0; offset42 = 0; offset43 = 0; offset44 = 0;



for k = 1:1:length(time)-1 
    
% phi1(k+1)
% phi21(k+1)
% phi31(k+1)
% phi41(k+1)
% phi1(k)
% phi31(k)

[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1), rdot_sim1(k+1), rdot_sim1_(k+1), phase1(k+1), phi1(k+1), phi21(k+1), phi31(k+1), phi41(k+1), offset11, offset21, offset31, offset41] = noisysim(x1,f1,Gt,M,X,PT, 0.7*GT,  GR,R,sigma,1,k,z,z_prev,phi1(k),phi31(k), time(k+1)-time(k), offset11, offset21, offset31, offset41);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1), rdot_sim2(k+1), rdot_sim2_(k+1), phase2(k+1), phi2(k+1), phi22(k+1), phi32(k+1), phi42(k+1), offset12, offset22, offset32, offset42] = noisysim(x2,f2,Gt,M,X,PT,   7*GT,  GR,R,sigma,2,k,z,z_prev,phi2(k),phi32(k), time(k+1)-time(k), offset12, offset22, offset32, offset42);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1), rdot_sim3(k+1), rdot_sim3_(k+1), phase3(k+1), phi3(k+1), phi23(k+1), phi33(k+1), phi43(k+1), offset13, offset23, offset33, offset43] = noisysim(x3,f3,Gt,M,X,PT,    GT,   GR,R,sigma,3,k,z,z_prev,phi3(k),phi33(k), time(k+1)-time(k), offset13, offset23, offset33, offset43);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1), rdot_sim4(k+1), rdot_sim4_(k+1), phase4(k+1), phi4(k+1), phi24(k+1), phi34(k+1), phi44(k+1), offset14, offset24, offset34, offset44] = noisysim(x4,f4,Gt,M,X,PT,0.5*GT,   GR,R,sigma,4,k,z,z_prev,phi4(k),phi34(k), time(k+1)-time(k), offset14, offset24, offset34, offset44);  
end
%%

figure
subplot(4,1,1),plot(time, r_sim1,'LineWidth',3);hold on; plot(time, r_sim1_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,2),plot(time, r_sim2,'LineWidth',3);hold on; plot(time, r_sim2_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,3),plot(time, r_sim3,'LineWidth',3);hold on; plot(time, r_sim3_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,4),plot(time, r_sim4,'LineWidth',3);hold on; plot(time, r_sim4_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
%%

figure
subplot(4,1,1),plot(time, rdot_sim1+1.6,'LineWidth',2);hold on; plot(time, radial(1,:),'LineWidth',2);title('Derived Radial Distance from Two Frequencies of Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,2),plot(time, rdot_sim2+1.9,'LineWidth',2);hold on; plot(time, radial(2,:),'LineWidth',2);title('Derived Radial Distance from Two Frequencies of Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,3),plot(time, 0.2*rdot_sim3+0.3,'LineWidth',2);hold on; plot(time, radial(3,:),'LineWidth',2);title('Derived Radial Distance from Two Frequencies of Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast')
subplot(4,1,4),plot(time, rdot_sim4+1,'LineWidth',2);hold on; plot(time, radial(4,:),'LineWidth',2);title('Derived Radial Distance from Two Frequencies of Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
%%
figure
subplot(4,1,1),plot(time, 1.6*rdot_sim1_+2.78,'LineWidth',3);hold on; plot(time, r_sim1_,'LineWidth',2);title('Derived Radial Distance from Two Frequencies Adjacent Time Stamps of Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,2),plot(time, 2*rdot_sim2_-0.03,'LineWidth',3);hold on; plot(time, r_sim2_,'LineWidth',2);title('Derived Radial Distance from Two Frequencies Adjacent Time Stamps of Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,3),plot(time, 2.4*rdot_sim3_+3.325,'LineWidth',3);hold on; plot(time, r_sim3_,'LineWidth',2);title('Derived Radial Distance from Two Frequencies Adjacent Time Stamps of Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')
subplot(4,1,4),plot(time, 1.8*rdot_sim4_+2.80,'LineWidth',3);hold on; plot(time, r_sim4_,'LineWidth',2);title('Derived Radial Distance from Two Frequencies Adjacent Time Stamps of Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth')

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

[magD12, magD22, magD32, magD42, phaseD12, phaseD22, phaseD32, phaseD42] = getMeas(time);

%%
figure 
subplot(411), 
plot(time, magD12, 'LineWidth', 2), title('Magnitude from Reader $\#1$ in 3D Motion','interpreter','latex'),hold on; plot(time, H1-0.000140, 'LineWidth', 2),
xlim([35, 155]);ylabel('Magnitude [V]'); legend('Measurement','Simulation','location','SouthEast'); grid on; grid minor


subplot(412), 
plot(time, magD22, 'LineWidth', 2), title('Magnitude from Reader $\#2$ in 3D Motion','interpreter','latex');hold on; plot(time, H2-0.00056, 'LineWidth', 2),
xlim([35, 155]);ylabel('Magnitude [V]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor

subplot(413), 
plot(time, magD32, 'LineWidth', 2), title('Magnitude from Reader $\#3$ in 3D Motion','interpreter','latex');hold on; plot(time, H3, 'LineWidth', 2),
xlim([35, 155]);ylabel('Magnitude [V]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor

subplot(414),
plot(time, magD42, 'LineWidth', 2), title('Magnitude from Reader $\#4$ in 3D Motion','interpreter','latex'),hold on; plot(time, H4-0.000110, 'LineWidth', 2),
xlim([35, 155]);legend('Measurement','Simulation','location','NorthEast'), ylabel('Magnitude [V]'); grid on; grid minor; xlabel('t [s]')

%%

figure
subplot(4,1,1),plot(time, phaseD12, 'LineWidth', 1),hold on; plot(time, 0.5*phase1+1.6,'LineWidth',2);hold on; %hold on; plot(time, r_sim1_,'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time, phaseD22, 'LineWidth', 1),hold on; plot(time, 0.5*phase2+1.6,'LineWidth',2);hold on; %hold on; plot(time, r_sim2_,'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time, phaseD32, 'LineWidth', 1),hold on; plot(time, 0.5*phase3+1.6,'LineWidth',2);hold on; %hold on; plot(time, r_sim3_,'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time, phaseD42, 'LineWidth', 1),hold on; plot(time, 0.5*phase4+1.6,'LineWidth',2);hold on; %hold on; plot(time, r_sim4_,'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

%%


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

