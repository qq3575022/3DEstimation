clc, clear, close all
%% =================================================================== Load Data ======================================================================
% ----------- Position of Four Readers ---------
x1 = [0,    0,    0.865];  
x2 = [2.29, 0,    1.27];   
x3 = [2.29, 2.52, 0.865]; 
x4 = [0,    2.52, 1.27];

% -------------- Time and Coordinates ----------
% time:     IMU Measurement
% coord3:   Ground Truth of 3D Coordinates
% z,z_prev: 3D coordinates for Simulation
[magD12, magD22, magD32, magD42, phaseD1, phaseD2, phaseD3, phaseD4, time] = getMeas();% Measurement Magnitude of Length 130854
%%

% -------------- Time and Coordinates ----------
% time:     IMU Measurement
% coord3:   Ground Truth of 3D Coordinates
% z,z_prev: 3D coordinates for Simulation

[coord3, radial, z, z_prev, H1, H2, H3, H4, H1_, H2_, H3_, H4_, r_sim, r_sim1, r_sim2, r_sim3, r_sim4, r_sim1_, r_sim2_, r_sim3_, r_sim4_, r_meas, r_meas1, r_meas2, r_meas3, r_meas4, r_phase, rphase1, rphase2, rphase3, rphase4, rdot_sim1, rdot_sim2, rdot_sim3, rdot_sim4, rdot_sim1_, rdot_sim2_, rdot_sim3_, rdot_sim4_, phase1, phase2, phase3, phase4, phigt1_1, phigt1_2, phigt1_3, phigt1_4, phigt2_1, phigt2_2, phigt2_3, phigt2_4, phi1_1, phi1_2, phi1_3, phi1_4, phi2_1, phi2_2, phi2_3, phi2_4, phi3_1, phi3_2, phi3_3, phi3_4, phi4_1, phi4_2, phi4_3, phi4_4] = get3Dcoord(x1, x2, x3, x4, time);
%%
% Start and End Index of 3D Motion
yST  = find(abs(time-107.99)<0.002);   yST = yST(1)-1;
yET  = find(abs(time-111.984)<0.002);  yET = yET(1)+1;

%
% Parameters of reader
Gt = 1/75*sqrt(1.462*3/4);     % tag's antenna gain
X  = 0.85;                     % polarization mismatch
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

% Parameters of reader
PT = 1;                             % reader's transmitted power
R = 15;
GT1 = 0.7*0.0331*sqrt(1.462*3/4);   % reader's trasmitter antenna gain -16.15dBi
GT2 =  7*0.0331*sqrt(1.462*3/4);    % reader's trasmitter antenna gain -6.15dBi
GT3 =    0.0331*sqrt(1.462*3/4);    % reader's trasmitter antenna gain -14.60dBi
GT4 = 0.5*0.0331*sqrt(1.462*3/4);   % reader's trasmitter antenna gain -17.61dBi

% Channel noise error covariance
sigma = 0.00012; 

% phase cconcatenation
offset11 = 0; offset12 = 0; offset13 = 0; offset14 = 0;  offset21 = 0; offset22 = 0; offset23 = 0; offset24 = 0;
offset31 = 0; offset32 = 0; offset33 = 0; offset34 = 0;  offset41 = 0; offset42 = 0; offset43 = 0; offset44 = 0;


for k = 1:1:length(time)-1 
[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1), rdot_sim1(k+1), rdot_sim1_(k+1), phase1(k+1), phigt1_1(k+1), phigt2_1(k+1), phi1_1(k+1), phi2_1(k+1), offset11, offset21, offset31, offset41] = noisysim(x1,f1,Gt,X,PT, GT1,R,sigma,0.005,k,z,z_prev,phigt1_1(k),phi1_1(k), time(k+1)-time(k), offset11, offset21, offset31, offset41);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1), rdot_sim2(k+1), rdot_sim2_(k+1), phase2(k+1), phigt1_2(k+1), phigt2_2(k+1), phi1_2(k+1), phi2_2(k+1), offset12, offset22, offset32, offset42] = noisysim(x2,f2,Gt,X,PT, GT2,R,sigma,0.0025,k,z,z_prev,phigt1_2(k),phi1_2(k), time(k+1)-time(k), offset12, offset22, offset32, offset42);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1), rdot_sim3(k+1), rdot_sim3_(k+1), phase3(k+1), phigt1_3(k+1), phigt2_3(k+1), phi1_3(k+1), phi2_3(k+1), offset13, offset23, offset33, offset43] = noisysim(x3,f3,Gt,X,PT, GT3, R,sigma,0.012,k,z,z_prev,phigt1_3(k),phi1_3(k), time(k+1)-time(k), offset13, offset23, offset33, offset43);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1), rdot_sim4(k+1), rdot_sim4_(k+1), phase4(k+1), phigt1_4(k+1), phigt2_4(k+1), phi1_4(k+1), phi2_4(k+1), offset14, offset24, offset34, offset44] = noisysim(x4,f4,Gt,X,PT, GT4, R,sigma,0.005,k,z,z_prev,phigt1_4(k),phi1_4(k), time(k+1)-time(k), offset14, offset24, offset34, offset44);  
end
%%
figure
subplot(4,1,1),plot(time, H1,'LineWidth',3);hold on; plot(time, H1_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);
subplot(4,1,2),plot(time, H2,'LineWidth',3);hold on; plot(time, H2_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);
subplot(4,1,3),plot(time, H3,'LineWidth',3);hold on; plot(time, H3_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth');ylim([1,2.5]),xlim([35, 155]);
subplot(4,1,4),plot(time, H4,'LineWidth',3);hold on; plot(time, H4_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);

%%
figure
subplot(4,1,1),plot(time, r_sim1,'LineWidth',3);hold on; plot(time, r_sim1_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);
subplot(4,1,2),plot(time, r_sim2,'LineWidth',3);hold on; plot(time, r_sim2_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);
subplot(4,1,3),plot(time, r_sim3,'LineWidth',3);hold on; plot(time, r_sim3_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth');ylim([1,2.5]),xlim([35, 155]);
subplot(4,1,4),plot(time, r_sim4,'LineWidth',3);hold on; plot(time, r_sim4_,'LineWidth',2);title('Derived Radial Distance from Magnitude of Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);
%%

figure
subplot(4,1,1),plot(time, rdot_sim1_,'LineWidth',3);hold on; plot(time, r_sim1_,'LineWidth',2);title('Derived Radial Distance from Two Frequencies of Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);
subplot(4,1,2),plot(time, rdot_sim2_,'LineWidth',3);hold on; plot(time, r_sim2_,'LineWidth',2);title('Derived Radial Distance from Two Frequencies of Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast'),xlim([35, 155]);
subplot(4,1,3),plot(time, rdot_sim3_,'LineWidth',3);hold on; plot(time, r_sim3_,'LineWidth',2);title('Derived Radial Distance from Two Frequencies of Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast'),xlim([35, 155]);
subplot(4,1,4),plot(time, rdot_sim4_,'LineWidth',3);hold on; plot(time, r_sim4_,'LineWidth',2);title('Derived Radial Distance from Two Frequencies of Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);

%%

figure
subplot(4,1,1),plot(time, rdot_sim1,'LineWidth',3);hold on; plot(time, radial(1,:),'LineWidth',2);title('Derived Radial Distance from Two Frequencies Adjacent Time Stamps of Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);ylim([1.5, 2.5]);
subplot(4,1,2),plot(time, rdot_sim2,'LineWidth',3);hold on; plot(time, radial(2,:),'LineWidth',2);title('Derived Radial Distance from Two Frequencies Adjacent Time Stamps of Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast'),xlim([35, 155]);ylim([1.2, 2.1]);
subplot(4,1,3),plot(time, rdot_sim3,'LineWidth',3);hold on; plot(time, radial(3,:),'LineWidth',2);title('Derived Radial Distance from Two Frequencies Adjacent Time Stamps of Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth','location','SouthEast'),xlim([35, 155]);ylim([0.8, 2.1]);
subplot(4,1,4),plot(time, rdot_sim4,'LineWidth',3);hold on; plot(time, radial(4,:),'LineWidth',2);title('Derived Radial Distance from Two Frequencies Adjacent Time Stamps of Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;legend('Simulation','Ground Truth'),xlim([35, 155]);ylim([1.1,2.1]);

%%
figure 
subplot(411), 
plot(time, magD12, 'LineWidth', 2), title('Magnitude from Reader $\#1$ in 3D Motion','interpreter','latex'),hold on; plot(time, H1, 'LineWidth', 2),%plot(time, H1-0.000140, 'LineWidth', 2),
xlim([35, 155]);ylabel('Magnitude [V]'); legend('Measurement','Simulation','location','SouthEast'); grid on; grid minor


subplot(412), 
plot(time, magD22, 'LineWidth', 2), title('Magnitude from Reader $\#2$ in 3D Motion','interpreter','latex');hold on; plot(time, H2, 'LineWidth', 2),%plot(time, H2-0.00056, 'LineWidth', 2),
xlim([35, 155]);ylabel('Magnitude [V]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor

subplot(413), 
plot(time, magD32, 'LineWidth', 2), title('Magnitude from Reader $\#3$ in 3D Motion','interpreter','latex');hold on; plot(time, H3, 'LineWidth', 2),
xlim([35, 155]);ylabel('Magnitude [V]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor,ylim([0, 1e-3])

subplot(414),
plot(time, magD42, 'LineWidth', 2), title('Magnitude from Reader $\#4$ in 3D Motion','interpreter','latex'),hold on; plot(time, H4, 'LineWidth', 2),%plot(time, H4-0.000110, 'LineWidth', 2),
xlim([35, 155]);legend('Measurement','Simulation','location','NorthEast'), ylabel('Magnitude [V]'); grid on; grid minor; xlabel('t [s]')

%%

figure
subplot(4,1,1),plot(time, phaseD1, 'LineWidth', 1),hold on; plot(time, 0.5*phase1+1.6,'LineWidth',2);hold on; 
title('Phase from Reader $\#1$ in 3D Motion','interpreter','latex'),xlim([35, 155]);ylabel('Phase [rad]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor%hold on; plot(time, r_sim1_,'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

subplot(4,1,2),plot(time, phaseD2, 'LineWidth', 1),hold on; plot(time, 0.5*phase2+1.6,'LineWidth',2);hold on; 
title('Phase from Reader $\#2$ in 3D Motion','interpreter','latex'),xlim([35, 155]);ylabel('Phase [rad]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor%hold on; plot(time, r_sim2_,'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

subplot(4,1,3),plot(time, phaseD3, 'LineWidth', 1),hold on; plot(time, 0.5*phase3+1.6,'LineWidth',2);hold on; 
title('Phase from Reader $\#3$ in 3D Motion','interpreter','latex'),xlim([35, 155]);ylabel('Phase [rad]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor%hold on; plot(time, r_sim3_,'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

subplot(4,1,4),plot(time, phaseD4, 'LineWidth', 1),hold on; plot(time, 0.5*phase4+1.6,'LineWidth',2);hold on; 
title('Phase from Reader $\#4$ in 3D Motion','interpreter','latex'),xlim([35, 155]);ylabel('Phase [rad]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor%hold on; plot(time, r_sim4_,'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
xlabel('t [s]')
%%
% 
% 
% [-mean(H1(2:124840) - magD12(2:124840)), rms(H1(2:124840) - magD12(2:124840)), sqrt(rms(H1(2:124840) - magD12(2:124840)))]
% [-mean(H2(2:124840) - magD22(2:124840)), rms(H2(2:124840) - magD22(2:124840)), sqrt(rms(H2(2:124840) - magD22(2:124840)))]
% [-mean(H3(2:124840) - magD32(2:124840)), rms(H3(2:124840) - magD32(2:124840)), sqrt(rms(H3(2:124840) - magD32(2:124840)))]
% [-mean(H4(2:124840) - magD42(2:124840)), rms(H4(2:124840) - magD42(2:124840)), sqrt(rms(H4(2:124840) - magD42(2:124840)))]
% 
% 
% [-mean(H1_(2:124840) - magD12(2:124840)), rms(H1_(2:124840) - magD12(2:124840)), sqrt(rms(H1_(2:124840) - magD12(2:124840)))]
% [-mean(H2_(2:124840) - magD22(2:124840)), rms(H2_(2:124840) - magD22(2:124840)), sqrt(rms(H2_(2:124840) - magD22(2:124840)))]
% [-mean(H3_(2:124840) - magD32(2:124840)), rms(H3_(2:124840) - magD32(2:124840)), sqrt(rms(H3_(2:124840) - magD32(2:124840)))]
% [-mean(H4_(2:124840) - magD42(2:124840)), rms(H4_(2:124840) - magD42(2:124840)), sqrt(rms(H4_(2:124840) - magD42(2:124840)))]
% 
% [-mean(H1_(2:124840) - H1(2:124840)), rms(H1_(2:124840) - H1(2:124840)), sqrt(rms(H1_(2:124840) - H1(2:124840)))]
% [-mean(H2_(2:124840) - H2(2:124840)), rms(H2_(2:124840) - H2(2:124840)), sqrt(rms(H2_(2:124840) - H2(2:124840)))]
% [-mean(H3_(2:124840) - H3(2:124840)), rms(H3_(2:124840) - H3(2:124840)), sqrt(rms(H3_(2:124840) - H3(2:124840)))]
% [-mean(H4_(2:124840) - H4(2:124840)), rms(H4_(2:124840) - H4(2:124840)), sqrt(rms(H4_(2:124840) - H4(2:124840)))]
% 
