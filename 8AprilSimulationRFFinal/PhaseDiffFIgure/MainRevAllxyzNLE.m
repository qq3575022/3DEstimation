clc, clear, close all

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
%Start and End Index of 3D Motion
yST  = find(abs(time-107.99)<0.002);   yST = yST(1)-1;
yET  = find(abs(time-111.984)<0.002);  yET = yET(1)+1;

%%
% Parameters of reader
Gt = 0.6*1.462*3/4;    % tag's antenna gain
X  = 0.85;             % polarization mismatch
M  = sqrt(4);          % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

% Parameters of reader
PT = 1;         % reader's transmitted power
GT = sqrt(1.462);     % reader's trasmitter antenna gain 9.5dBi
GR = sqrt(1.462);     % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.0012; 

offset11 = 0; offset21 = 0; offset31 = 0; offset41 = 0;offset12 = 0; offset22 = 0; offset32 = 0; offset42 = 0;
offset13 = 0; offset23 = 0; offset33 = 0; offset43 = 0;offset14 = 0; offset24 = 0; offset34 = 0; offset44 = 0;
offset51 = 0; offset52 = 0; offset53 = 0; offset54 = 0;offset61 = 0; offset62 = 0; offset63 = 0; offset64 = 0; 


for k = yST-1:1:yET  

[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1),r_meas1(k+1),rphase1(k+1),rdot_sim1(k+1),rdot_sim1_(k+1),phigt1_1(k+1),phigt2_1(k+1),phi1_1(k+1),phi2_1(k+1),offset11,offset21,offset31,offset41,offset51,offset61] = noisysimNoMultiStatic(x1,f1,Gt,M,X,PT,0.7*GT,GR,R,sigma,k,z,z_prev,phigt1_1(k),phi1_1(k),time(k)-time(k-1),magD12(k),offset11,offset21,offset31,offset41,offset51,offset61);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1),r_meas2(k+1),rphase2(k+1),rdot_sim2(k+1),rdot_sim2_(k+1),phigt1_2(k+1),phigt2_2(k+1),phi1_2(k+1),phi2_2(k+1),offset12,offset22,offset32,offset42,offset52,offset62] = noisysimNoMultiStatic(x2,f2,Gt,M,X,PT,  7*GT,GR,R,sigma,k,z,z_prev,phigt1_2(k),phi1_2(k),time(k)-time(k-1),magD22(k),offset12,offset22,offset32,offset42,offset52,offset62);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1),r_meas3(k+1),rphase3(k+1),rdot_sim3(k+1),rdot_sim3_(k+1),phigt1_3(k+1),phigt2_3(k+1),phi1_3(k+1),phi2_3(k+1),offset13,offset23,offset33,offset43,offset53,offset63] = noisysimNoMultiStatic(x3,f3,Gt,M,X,PT,    GT,GR,R,sigma,k,z,z_prev,phigt1_3(k),phi1_3(k),time(k)-time(k-1),magD32(k),offset13,offset23,offset33,offset43,offset53,offset63);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1),r_meas4(k+1),rphase4(k+1),rdot_sim4(k+1),rdot_sim4_(k+1),phigt1_4(k+1),phigt2_4(k+1),phi1_4(k+1),phi2_4(k+1),offset14,offset24,offset34,offset44,offset54,offset64] = noisysimNoMultiStatic(x4,f4,Gt,M,X,PT,0.5*GT,GR,R,sigma,k,z,z_prev,phigt1_4(k),phi1_4(k),time(k)-time(k-1),magD42(k),offset14,offset24,offset34,offset44,offset54,offset64);  

end
%%
figure
subplot(4,1,1),plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

%%
figure
subplot(4,1,1),plot(time(yST:yET), r_meas1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), r_meas2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), r_meas3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), r_meas4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
%%
figure
subplot(4,1,1),plot(time(yST:yET), r_sim1_(yST:yET),'LineWidth',3);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',1);title('Simulated Radial Distance Derived from Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([1, 3]);legend('Simulation','Ground truth');
subplot(4,1,2),plot(time(yST:yET), r_sim2_(yST:yET),'LineWidth',3);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',1);title('Simulated Radial Distance Derived from Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([0.8, 2.8]);legend('Simulation','Ground truth');
subplot(4,1,3),plot(time(yST:yET), r_sim3_(yST:yET),'LineWidth',3);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',1);title('Simulated Radial Distance Derived from Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([0.5, 2.5]);legend('Simulation','Ground truth');
subplot(4,1,4),plot(time(yST:yET), r_sim4_(yST:yET),'LineWidth',3);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',1);title('Simulated Radial Distance Derived from Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([0.8, 2.8]);legend('Simulation','Ground truth');
%%
figure
subplot(4,1,1),plot(time(yST:yET), rphase1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,2),plot(time(yST:yET), rphase2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,3),plot(time(yST:yET), rphase3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,4),plot(time(yST:yET), rphase4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');


%%

%%

% measrdot1 = phaseDiff1; measrdot2 = phaseDiff2; measrdot3 = phaseDiff3; measrdot4 = phaseDiff4; 
% 
% dmeasr1 = zeros(1, length(time));dmeasr2 = zeros(1, length(time));
% dmeasr3 = zeros(1, length(time));dmeasr4 = zeros(1, length(time));
% 
% dmeasrdot1 = zeros(1, length(time));  dmeasrdot2 = zeros(1, length(time)); 
% dmeasrdot3 = zeros(1, length(time));  dmeasrdot4 = zeros(1, length(time)); 
% 
% for i = 3:1:length(time)
%     difT = time(i) - time(1);
%     index = min(round(difT*(4.8469*1000)), length(measrdot1));
%     index2 = min(round(difT*(9.8607*1000)), length(measrdot3));
%     dmeasr1(i) = measr1(index);
%     dmeasr2(i) = measr2(index);
%     dmeasr3(i) = measr3(index2);
%     dmeasr4(i) = measr4(index);
%     
%     dmeasrdot1(i) = measrdot1(index);
%     dmeasrdot2(i) = measrdot2(index);
%     dmeasrdot3(i) = measrdot3(index2);
%     dmeasrdot4(i) = measrdot4(index);
% end
% 
% dmeasr1(1) = measr1(1);dmeasr2(1) = measr2(1);dmeasr3(1) = measr3(1);dmeasr4(1) = measr4(1);
% dmeasrdot1(1) = measrdot1(1);dmeasrdot2(1) = measrdot2(1);dmeasrdot3(1) = measrdot3(1);dmeasrdot4(1) = measrdot4(1);
% dmeasr1(2) = measr1(2);dmeasr2(2) = measr2(2);dmeasr3(2) = measr3(2);dmeasr4(2) = measr4(2);
% dmeasrdot1(2) = measrdot1(2);dmeasrdot2(2) = measrdot2(2);dmeasrdot3(2) = measrdot3(2);dmeasrdot4(2) = measrdot4(2);

%%
figure
subplot(411), plot(time,  r_meas1, 'LineWidth', 1), hold on, plot(time, r_sim1, 'LineWidth', 1), title('Radial Distance $r1$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);
subplot(412), plot(time,  r_meas2, 'LineWidth', 1), hold on, plot(time, r_sim2, 'LineWidth', 1),title('Radial Distance $r2$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);
subplot(413), plot(time,  r_meas3, 'LineWidth', 1), hold on, plot(time, r_sim3, 'LineWidth', 1),title('Radial Distance $r3$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);
subplot(414), plot(time,  r_meas4, 'LineWidth', 1), hold on, plot(time, r_sim4, 'LineWidth', 1),title('Radial Distance $r4$','interpreter','latex'); xlabel('time[s]'), ylabel('Distance [m]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);

%%
figure
subplot(411), plot(time(yST:yET), phaseD1(yST+1300:yET+1300), 'LineWidth', 2); hold on; plot(time(yST:yET), rdot_sim1(yST:yET), 'LineWidth', 2), title('Radial Velocity $\dot r_1$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);grid on; grid minor;
subplot(412), plot(time(yST:yET), phaseD2(yST+1300:yET+1300), 'LineWidth', 2); hold on; plot(time(yST:yET), rdot_sim2(yST:yET), 'LineWidth', 2), title('Radial Velocity $\dot r_2$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);grid on; grid minor;
subplot(413), plot(time(yST:yET), phaseD3(yST+1300:yET+1300), 'LineWidth', 2); hold on; plot(time(yST:yET), rdot_sim3(yST:yET), 'LineWidth', 2), title('Radial Velocity $\dot r_3$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth','Location','SouthEast'), xlim([107.99, 111.984]);grid on; grid minor;
subplot(414), plot(time(yST:yET), phaseD4(yST+1300:yET+1300), 'LineWidth', 2); hold on; plot(time(yST:yET), rdot_sim4(yST:yET), 'LineWidth', 2), title('Radial Velocity $\dot r_4$','interpreter','latex'); xlabel('time[s]'), ylabel('Velocity [m/s]','interpreter','latex'), legend('measurement','ground truth'), xlim([107.99, 111.984]);grid on; grid minor;

