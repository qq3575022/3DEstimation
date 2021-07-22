clc, clear, close all
% =================================================================== Load Data ======================================================================
% ----------- Position of Four Readers ---------
x1 = [0,    0,    0.865];  
x2 = [2.29, 0,    1.27];   
x3 = [2.29, 2.52, 0.865]; 
x4 = [0,    2.52, 1.27];

% -------------- Time and Coordinates ----------
% time:     IMU Measurement
% coord3:   Ground Truth of 3D Coordinates
% z,z_prev: 3D coordinates for Simulation
[magD12, magD22, magD32, magD42, time] = getMeas();% Measurement Magnitude of Length 130854
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
%%
% =================================================================== Simulation ======================================================================

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
sigma = 0.02; 

offset11 = 0; offset21 = 0; offset31 = 0; offset41 = 0;offset12 = 0; offset22 = 0; offset32 = 0; offset42 = 0;
offset13 = 0; offset23 = 0; offset33 = 0; offset43 = 0;offset14 = 0; offset24 = 0; offset34 = 0; offset44 = 0;
offset51 = 0; offset52 = 0; offset53 = 0; offset54 = 0;offset61 = 0; offset62 = 0; offset63 = 0; offset64 = 0; 

Kfac = 0.05;
rng(2,'twister');

for k = yST-1:1:yET  

[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1),r_meas1(k+1),rphase1(k+1),rdot_sim1(k+1),rdot_sim1_(k+1),phigt1_1(k+1),phigt2_1(k+1),phi1_1(k+1),phi2_1(k+1),offset11,offset21,offset31,offset41,offset51,offset61] = noisysimNoMultiStatic(x1,f1,Gt,X,PT,GT1,R,sigma,0.025,k,z,z_prev,phigt1_1(k),phi1_1(k),time(k)-time(k-1),magD12(k),offset11,offset21,offset31,offset41,offset51,offset61);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1),r_meas2(k+1),rphase2(k+1),rdot_sim2(k+1),rdot_sim2_(k+1),phigt1_2(k+1),phigt2_2(k+1),phi1_2(k+1),phi2_2(k+1),offset12,offset22,offset32,offset42,offset52,offset62] = noisysimNoMultiStatic(x2,f2,Gt,X,PT,GT2,R,sigma,0.0025,k,z,z_prev,phigt1_2(k),phi1_2(k),time(k)-time(k-1),magD22(k),offset12,offset22,offset32,offset42,offset52,offset62);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1),r_meas3(k+1),rphase3(k+1),rdot_sim3(k+1),rdot_sim3_(k+1),phigt1_3(k+1),phigt2_3(k+1),phi1_3(k+1),phi2_3(k+1),offset13,offset23,offset33,offset43,offset53,offset63] = noisysimNoMultiStatic(x3,f3,Gt,X,PT,GT3,R,sigma,0.088,k,z,z_prev,phigt1_3(k),phi1_3(k),time(k)-time(k-1),magD32(k),offset13,offset23,offset33,offset43,offset53,offset63);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1),r_meas4(k+1),rphase4(k+1),rdot_sim4(k+1),rdot_sim4_(k+1),phigt1_4(k+1),phigt2_4(k+1),phi1_4(k+1),phi2_4(k+1),offset14,offset24,offset34,offset44,offset54,offset64] = noisysimNoMultiStatic(x4,f4,Gt,X,PT,GT4,R,sigma,0.018,k,z,z_prev,phigt1_4(k),phi1_4(k),time(k)-time(k-1),magD42(k),offset14,offset24,offset34,offset44,offset54,offset64);  


rphase1(k+1) = rphase1(k+1)+0.0;
rphase2(k+1) = rphase2(k+1)+0.0;
rphase3(k+1) = rphase3(k+1)+0.0;
rphase4(k+1) = rphase4(k+1)+0.0;

[r_sim(1,k+1),r_sim(2,k+1),r_sim(3,k+1)]   = trilateration(x1, x2, x3, x4, r_sim1(k+1), r_sim2(k+1), r_sim3(k+1),  r_sim4(k+1));
[r_phase(1,k+1),r_phase(2,k+1),r_phase(3,k+1)] = trilateration(x1, x2, x3, x4, rphase1(k+1),rphase2(k+1),rphase3(k+1),rphase4(k+1));
[r_meas(1,k+1),r_meas(2,k+1),r_meas(3,k+1)]  = trilateration(x1, x2, x3, x4, r_meas1(k+1),r_meas2(k+1),r_meas3(k+1),r_meas4(k+1));

end

rtime = time(yST:yET); 

r1 = rphase1(yST:yET); r2 = rphase2(yST:yET); r3 = rphase3(yST:yET); r4 = rphase4(yST:yET);
rdot1 = rdot_sim1(yST:yET); rdot2 = rdot_sim2(yST:yET); rdot3 = rdot_sim3(yST:yET); rdot4 = rdot_sim4(yST:yET);

%save('SimRFPhase1.mat', 'r1', 'r2', 'r3', 'r4', 'rdot1', 'rdot2', 'rdot3', 'rdot4','rtime');
r1 = r_sim1(yST:yET); r2 = r_sim2(yST:yET); r3 = r_sim3(yST:yET); r4 = r_sim4(yST:yET);

%save('SimRFMag1.mat', 'r1', 'r2', 'r3', 'r4', 'rdot1', 'rdot2', 'rdot3', 'rdot4','rtime');


%% ==================================================================== Results ======================================================================
% -------------- Radial Distances RMS: Reader1 Reader2 Reader3 Reader4 ----------
fprintf('================================== radial distance RMS ================================== ')
rmeas  = 1000*[rms(r_meas1(yST:yET) - radial(1,yST:yET)), rms(r_meas2(yST:yET) - radial(2,yST:yET)), rms(r_meas3(yST:yET) - radial(3,yST:yET)), rms(r_meas4(yST:yET) - radial(4,yST:yET))]
rmag   = 1000*[rms(r_sim1(yST:yET)  - radial(1,yST:yET)), rms(r_sim2(yST:yET)  - radial(2,yST:yET)), rms(r_sim3(yST:yET) - radial(3,yST:yET)),  rms(r_sim4(yST:yET)  - radial(4,yST:yET))]
rphase = 1000*[rms(rphase1(yST:yET) - radial(1,yST:yET)), rms(rphase2(yST:yET) - radial(2,yST:yET)), rms(rphase3(yST:yET) - radial(3,yST:yET)), rms(rphase4(yST:yET) - radial(4,yST:yET))]

fprintf('================================== x y z triangulation RMS in total ================================== ')
meas  = 1000*sqrt(rms(coord3(1,yST:yET) - r_meas(1,yST:yET)).^2  + rms(coord3(2,yST:yET) - r_meas(2,yST:yET)).^2  + rms(coord3(3,yST:yET) - r_meas(3,yST:yET)).^2)
mag   = 1000*sqrt(rms(coord3(1,yST:yET) - r_sim(1,yST:yET)).^2   + rms(coord3(2,yST:yET) - r_sim(2,yST:yET)).^2   + rms(coord3(3,yST:yET) - r_sim(3,yST:yET)).^2)
phase = 1000*sqrt(rms(coord3(1,yST:yET) - r_phase(1,yST:yET)).^2 + rms(coord3(2,yST:yET) - r_phase(2,yST:yET)).^2 + rms(coord3(3,yST:yET) - r_phase(3,yST:yET)).^2)

fprintf('================================== x y triangulation RMS in total ================================== ')
meas  = 1000*sqrt(rms(coord3(1,yST:yET) - r_meas(1,yST:yET)).^2  + rms(coord3(2,yST:yET) - r_meas(2,yST:yET)).^2)
mag   = 1000*sqrt(rms(coord3(1,yST:yET) - r_sim(1,yST:yET)).^2   + rms(coord3(2,yST:yET) - r_sim(2,yST:yET)).^2)
phase = 1000*sqrt(rms(coord3(1,yST:yET) - r_phase(1,yST:yET)).^2 + rms(coord3(2,yST:yET) - r_phase(2,yST:yET)).^2 )

%%
fprintf('================================== x y z triangulation RMS separately ================================== ')
meas  = 1000*[rms(coord3(1,yST:yET) - r_meas(1,yST:yET)), rms(coord3(2,yST:yET) - r_meas(2,yST:yET)),  rms(coord3(3,yST:yET) - r_meas(3,yST:yET))]
mag   = 1000*[rms(coord3(1,yST:yET) - r_sim(1,yST:yET)),  rms(coord3(2,yST:yET) - r_sim(2,yST:yET)),   rms(coord3(3,yST:yET) - r_sim(3,yST:yET))]
phase = 1000*[rms(coord3(1,yST:yET) - r_phase(1,yST:yET)),rms(coord3(2,yST:yET) - r_phase(2,yST:yET)), rms(coord3(3,yST:yET) - r_phase(3,yST:yET))]

%%
% fprintf('================================== radial distance mean value ================================== ')
% rmeas  = 100*[mean(r_meas1(yST:yET) - radial(1,yST:yET)), mean(r_meas2(yST:yET) - radial(2,yST:yET)), mean(r_meas3(yST:yET) - radial(3,yST:yET)), mean(r_meas4(yST:yET) - radial(4,yST:yET))]
% rmag   = 100*[mean(r_sim1(yST:yET)  - radial(1,yST:yET)), mean(r_sim2(yST:yET)  - radial(2,yST:yET)), mean(r_sim3(yST:yET) - radial(3,yST:yET)),  mean(r_sim4(yST:yET)  - radial(4,yST:yET))]
% rphase = 100*[mean(rphase1(yST:yET) - radial(1,yST:yET)), mean(rphase2(yST:yET) - radial(2,yST:yET)), mean(rphase3(yST:yET) - radial(3,yST:yET)), mean(rphase4(yST:yET) - radial(4,yST:yET))]
% 
% fprintf('================================== x y z triangulation Mean Error in total ================================== ')
% meas  = 100*sqrt(mean(coord3(1,yST:yET) - r_meas(1,yST:yET)).^2  + mean(coord3(2,yST:yET) - r_meas(2,yST:yET)).^2  + mean(coord3(3,yST:yET) - r_meas(3,yST:yET)).^2)
% mag   = 100*sqrt(mean(coord3(1,yST:yET) - r_sim(1,yST:yET)).^2   + mean(coord3(2,yST:yET) - r_sim(2,yST:yET)).^2   + mean(coord3(3,yST:yET) - r_sim(3,yST:yET)).^2)
% phase = 100*sqrt(mean(coord3(1,yST:yET) - r_phase(1,yST:yET)).^2 + mean(coord3(2,yST:yET) - r_phase(2,yST:yET)).^2 + mean(coord3(3,yST:yET) - r_phase(3,yST:yET)).^2)

%%
% -------------- Measurement Radial Distances: Reader1 Reader2 Reader3 Reader4 ----------

figure
subplot(4,1,1),plot(time(yST:yET), r_meas1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim1__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time(yST:yET), r_meas2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim2__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time(yST:yET), r_meas3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim3__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time(yST:yET), r_meas4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);%hold on; plot(time(yST:yET), r_sim4__(yST:yET),'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

% -------------- Magnitude Radial Distances: Reader1 Reader2 Reader3 Reader4 ----------
%%
figure
subplot(4,1,1),plot(time(yST:yET), r_sim1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([1, 3]);legend('Simulation','Ground truth');
subplot(4,1,2),plot(time(yST:yET), r_sim2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([0.8, 2.8]);legend('Simulation','Ground truth');
subplot(4,1,3),plot(time(yST:yET), r_sim3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([0.5, 2.5]);legend('Simulation','Ground truth');
subplot(4,1,4),plot(time(yST:yET), r_sim4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);ylim([0.8, 2.8]);legend('Simulation','Ground truth');
%%
figure
subplot(4,1,1),plot(time(yST:yET), rdot_sim1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), rdot_sim1_(yST:yET),'LineWidth',2);title('Simulated Radial Velocity Derived from Phase Difference from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Velocity [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');%ylim([1, 3]);
subplot(4,1,2),plot(time(yST:yET), rdot_sim2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), rdot_sim2_(yST:yET),'LineWidth',2);title('Simulated Radial Velocity Derived from Phase Difference from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Velocity [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');%ylim([0.8, 2.8]);legend('Simulation','Ground truth');
subplot(4,1,3),plot(time(yST:yET), rdot_sim3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), rdot_sim3_(yST:yET),'LineWidth',2);title('Simulated Radial Velocity Derived from Phase Difference from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Velocity [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');%ylim([0.5, 2.5]);legend('Simulation','Ground truth');
subplot(4,1,4),plot(time(yST:yET), rdot_sim4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), rdot_sim4_(yST:yET),'LineWidth',2);title('Simulated Radial Velocity Derived from Phase Difference from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Velocity [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');%ylim([0.8, 2.8]);legend('Simulation','Ground truth');

% -------------- Two Frequencies Radial Distances: Reader1 Reader2 Reader3 Reader4 ----------
%%
figure
subplot(4,1,1),plot(time(yST:yET), rphase1(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(1,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,2),plot(time(yST:yET), rphase2(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(2,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,3),plot(time(yST:yET), rphase3(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(3,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(4,1,4),plot(time(yST:yET), rphase4(yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), radial(4,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#4$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');

%%
figure
subplot(3,1,1),plot(time(yST:yET), r_phase(1,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(1,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(3,1,2),plot(time(yST:yET), r_phase(2,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(2,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(3,1,3),plot(time(yST:yET), r_phase(3,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(3,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');

%%
figure
subplot(3,1,1),plot(time(yST:yET), r_sim(1,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(1,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#1$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(3,1,2),plot(time(yST:yET), r_sim(2,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(2,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#2$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
subplot(3,1,3),plot(time(yST:yET), r_sim(3,yST:yET),'LineWidth',2);hold on; plot(time(yST:yET), coord3(3,yST:yET),'LineWidth',2);title('Simulated Radial Distance Derived from Phase from Reader $\#3$ in 3D','interpreter','latex');ylabel('Radial Distance [m]');xlabel('t [s]');grid on; grid minor;xlim([108, 112]);legend('Simulation','Ground truth');
