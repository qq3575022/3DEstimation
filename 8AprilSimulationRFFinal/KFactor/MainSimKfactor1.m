clc, clear, close all

% ================================ Get  coordinates ==============================

[time, coord3, z, z_prev, H1, H2, H3, H4, H1_, H2_, H3_, H4_, r_sim1, r_sim2, r_sim3, r_sim4, r_sim1_, r_sim2_, r_sim3_, r_sim4_, rdot_sim1, rdot_sim2, rdot_sim3, rdot_sim4, rdot_sim1_, rdot_sim2_, rdot_sim3_, rdot_sim4_, diff1, diff2, diff3, diff4] = get3Dcoord();


% ============================== Simulation parameter ==============================

Gt = 1/75*sqrt(1.462*3/4);    % tag's antenna gain
X = 0.85;      % polarization mismatch
M = sqrt(4);         % load modulation factor of the tag
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
sigma = 0.00012; 

% ================================= Reader location =================================

x1 = [0,        0,   0.865];  
x2 = [2.29,     0,    1.27];   
x3 = [2.29,  2.52,   0.865]; 
x4 = [0,     2.52,    1.27];

% phase cconcatenation
% global l1; global l2; global l3; l1 = 0;l2 = 0;l3 = 0;k = 1; 
multi = randi(10, [20,1]);

for k = 1:1:length(time)-1  
    
[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1), rdot_sim1(k+1), rdot_sim1_(k+1)] = noisysim5(x1,f1,Gt,M,X,PT,0.7*GT,GR, R,sigma,1,k,z,z_prev,time(k+1)-time(k),time(k),multi);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1), rdot_sim2(k+1), rdot_sim2_(k+1)] = noisysim5(x2,f2,Gt,M,X,PT,7*GT,GR, R,sigma,2,k,z,z_prev,time(k+1)-time(k),time(k),multi);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1), rdot_sim3(k+1), rdot_sim3_(k+1)] = noisysim5(x3,f3,Gt,M,X,PT,GT,GR, R,sigma,3,k,z,z_prev,time(k+1)-time(k),time(k),multi);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1), rdot_sim4(k+1), rdot_sim4_(k+1)] = noisysim5(x4,f4,Gt,M,X,PT,0.5*GT,GR, R,sigma,4,k,z,z_prev,time(k+1)-time(k),time(k),multi);  

end

r_sim_gt = NaN(4,length(coord3));

for i = 1:1:length(coord3)
    r_sim_gt(1,i) = sqrt((coord3(1,i) - x1(1))^2 + (coord3(2,i) - x1(2))^2 +(coord3(3,i) - x1(3))^2);
    r_sim_gt(2,i) = sqrt((coord3(1,i) - x2(1))^2 + (coord3(2,i) - x2(2))^2 +(coord3(3,i) - x2(3))^2);
    r_sim_gt(3,i) = sqrt((coord3(1,i) - x3(1))^2 + (coord3(2,i) - x3(2))^2 +(coord3(3,i) - x3(3))^2);
    r_sim_gt(4,i) = sqrt((coord3(1,i) - x4(1))^2 + (coord3(2,i) - x4(2))^2 +(coord3(3,i) - x4(3))^2);

    
end
%%
figure
subplot(4,1,1),plot(time, H1,'LineWidth',2);hold on; plot(time, H1_,'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time, H2,'LineWidth',2);hold on; plot(time, H2_,'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time, H3,'LineWidth',2);hold on; plot(time, H3_,'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time, H4,'LineWidth',2);hold on; plot(time, H4_,'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

figure
subplot(4,1,1),plot(time, r_sim1,'LineWidth',2);hold on; plot(time, r_sim_gt(1,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#1$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time, r_sim2,'LineWidth',2);hold on; plot(time, r_sim_gt(2,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#2$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time, r_sim3,'LineWidth',2);hold on; plot(time, r_sim_gt(3,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#3$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time, r_sim4,'LineWidth',2);hold on; plot(time, r_sim_gt(4,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#4$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
%%
figure
subplot(4,1,1),plot(time, rdot_sim1_,'LineWidth',2);hold on; plot(time, rdot_sim1_,'LineWidth',2);title('Simulated Radial Distance from Reader $\#1$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time, rdot_sim2_,'LineWidth',2);hold on; plot(time, rdot_sim2_,'LineWidth',2);title('Simulated Radial Distance from Reader $\#2$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time, rdot_sim3_,'LineWidth',2);hold on; plot(time, rdot_sim3_,'LineWidth',2);title('Simulated Radial Distance from Reader $\#3$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time, rdot_sim4_,'LineWidth',2);hold on; plot(time, rdot_sim4_,'LineWidth',2);title('Simulated Radial Distance from Reader $\#4$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;

%%
%%
figure
subplot(4,1,1),plot(time, rdot_sim1,'LineWidth',2);hold on; plot(time, r_sim_gt(1,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#1$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time, rdot_sim2,'LineWidth',2);hold on; plot(time, r_sim_gt(2,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#2$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time, rdot_sim3,'LineWidth',2);hold on; plot(time, r_sim_gt(3,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#3$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time, rdot_sim4,'LineWidth',2);hold on; plot(time, r_sim_gt(4,:),'LineWidth',2);title('Simulated Radial Distance from Reader $\#4$ in 3D','interpreter','latex');legend('sim','gt');xlabel('t [s]');grid on; grid minor;

%%
%
%
[magD12, magD22, magD32, magD42] = getMeas(time);
%%
figure 
subplot(411), 
plot(time, magD12, 'LineWidth', 2), title('Reader 1 Magnitude'),hold on; plot(time, H1-0.000140, 'LineWidth', 2),
xlim([35, 155]);ylabel('Magnitude [V]'); legend('Measurement','Simulation','location','SouthEast'); grid on; grid minor


subplot(412), 
plot(time, magD22, 'LineWidth', 2), title('Reader 2 Magnitude');hold on; plot(time, H2-0.00056, 'LineWidth', 2),
xlim([35, 155]);ylabel('Magnitude [V]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor

subplot(413), 
plot(time, magD32, 'LineWidth', 2), title('Reader 3 Magnitude');hold on; plot(time, H3, 'LineWidth', 2),
xlim([35, 155]);ylabel('Magnitude [V]'); legend('Measurement','Simulation','location','NorthEast'); grid on; grid minor

subplot(414),
plot(time, magD42, 'LineWidth', 2), title('Reader 4 Magnitude'),hold on; plot(time, H4-0.000110, 'LineWidth', 2),
xlim([35, 155]);legend('Measurement','Simulation','location','NorthEast'), ylabel('Magnitude [V]'); grid on; grid minor



[-mean(H1(2:124840) - magD12(2:124840)),  rms(H1(2:124840) - magD12(2:124840)),  sqrt(rms(H1(2:124840) - magD12(2:124840)))]
[-mean(H2(2:124840) - magD22(2:124840)),  rms(H2(2:124840) - magD22(2:124840)),  sqrt(rms(H2(2:124840) - magD22(2:124840)))]
[-mean(H3(2:124840) - magD32(2:124840)),  rms(H3(2:124840) - magD32(2:124840)),  sqrt(rms(H3(2:124840) - magD32(2:124840)))]
[-mean(H4(2:124840) - magD42(2:124840)),  rms(H4(2:124840) - magD42(2:124840)),  sqrt(rms(H4(2:124840) - magD42(2:124840)))]


[-mean(H1_(2:124840) - magD12(2:124840)), rms(H1_(2:124840) - magD12(2:124840)), sqrt(rms(H1_(2:124840) - magD12(2:124840)))]
[-mean(H2_(2:124840) - magD22(2:124840)), rms(H2_(2:124840) - magD22(2:124840)), sqrt(rms(H2_(2:124840) - magD22(2:124840)))]
[-mean(H3_(2:124840) - magD32(2:124840)), rms(H3_(2:124840) - magD32(2:124840)), sqrt(rms(H3_(2:124840) - magD32(2:124840)))]
[-mean(H4_(2:124840) - magD42(2:124840)), rms(H4_(2:124840) - magD42(2:124840)), sqrt(rms(H4_(2:124840) - magD42(2:124840)))]

[-mean(H1_(2:124840) - H1(2:124840)),     rms(H1_(2:124840) - H1(2:124840)),     sqrt(rms(H1_(2:124840) - H1(2:124840)))]
[-mean(H2_(2:124840) - H2(2:124840)),     rms(H2_(2:124840) - H2(2:124840)),     sqrt(rms(H2_(2:124840) - H2(2:124840)))]
[-mean(H3_(2:124840) - H3(2:124840)),     rms(H3_(2:124840) - H3(2:124840)),     sqrt(rms(H3_(2:124840) - H3(2:124840)))]
[-mean(H4_(2:124840) - H4(2:124840)),     rms(H4_(2:124840) - H4(2:124840)),     sqrt(rms(H4_(2:124840) - H4(2:124840)))]
%%

%Get K Value
m41 = getm3(sqrt(magD12));
m42 = getm3(sqrt(magD22));
m43 = getm3(sqrt(magD32));
m44 = getm3(sqrt(magD42));
%%
figure
plot(time(1:end-10000), medfilt1(m41, 10000),'LineWidth', 3);title('K Value for Reader 1 2 3 4 in Measurement');ylabel('K Value'); grid on; grid minor;xlim([5,130]);ylim([2,16000]);%ylim([0,1000]);xlim([5,140])
hold on, plot(time(1:end-10000), medfilt1(m42, 10000),'LineWidth', 3);%title('K Value for Reader 2 in Simulation');ylabel('K Value'); grid on; grid minor;%ylim([0,1000]);xlim([5,140])
hold on, plot(time(1:end-10000), medfilt1(m43, 10000),'LineWidth', 3);%title('K Value for Reader 3 in Simulation');ylabel('K Value'); grid on; grid minor;%ylim([0,1000]);xlim([5,140])
hold on, plot(time(1:end-10000), medfilt1(m44, 10000),'LineWidth', 3);legend('Reader 1','Reader 2', 'Reader 3', 'Reader 4');xlabel('time [s]');%title('K Value for Reader 4 in Simulation');ylabel('K Value'); grid on; grid minor;xlabel('time [s]');%ylim([0,1000]);xlim([5,140]);xlabel('time [s]');
%set(gca, 'YScale', 'log')



% ----------------------------------------------- End ========================================






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

