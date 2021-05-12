clc, clear, close all

[time, coord3] = get3Dcoord();

%
Gt = 1.462*3/4;    % tag's antenna gain
X = 0.85;          % polarization mismatch
M = 4;             % load modulation factor of the tag
f1 = 5.8*10^9;
f2 = 5.83*10^9;
f3 = 5.82*10^9;
f4 = 5.85*10^9;

x1 = [0,    0,  0.865];  %[2.6256, 0.0889,0.858];%[-0.2,-0.2, 0.4];
x2 = [2.29, 0,  1.27];   %[ 0.9, 0.8, 0.2];
x3 = [2.29,2.52, 0.865]; %[ 0.8,-0.3,-0.2];
x4 = [0, 2.52,  1.27];

% Parameters of reader
PT = 1;         % reader's transmitted power
GT = 1.462;     % reader's trasmitter antenna gain 9.5dBi
GR = 1.462;     % reader's receiver   antenna gain 9.5dBi
R = 15;

% Channel noise error covariance
sigma = 0.0000012; 

% phase cconcatenation
% global l1; global l2; global l3; l1 = 0;l2 = 0;l3 = 0;k = 1; 

% Get RSS and phase from each reader observing the moving tag
z = NaN(3,length(coord3)-1); z_prev = NaN(3,length(coord3)-1);

z_prev(1,:) = coord3(1,1:end-1); z(1,:) = coord3(1,2:end);% x coordinate
z_prev(2,:) = coord3(2,1:end-1); z(2,:) = coord3(2,2:end);% y coordinate
z_prev(3,:) = coord3(3,1:end-1); z(3,:) = coord3(3,2:end);% z coordinate

zdiff1 = z(1,:)-z_prev(1,:);
zdiff2 = z(2,:)-z_prev(2,:);
zdiff3 = z(3,:)-z_prev(3,:);

%
H1 = NaN(1,length(coord3));   H1_ = NaN(1,length(coord3));  %phi1 = NaN(1,length(coord3));       phi_mu1 = NaN(1,length(coord3));
H2 = NaN(1,length(coord3));   H2_ = NaN(1,length(coord3));   %phi2 = NaN(1,length(coord3));       phi_mu2 = NaN(1,length(coord3));
H3 = NaN(1,length(coord3));   H3_ = NaN(1,length(coord3)); %phi3 = NaN(1,length(coord3));       phi_mu3 = NaN(1,length(coord3));
H4 = NaN(1,length(coord3));   H4_ = NaN(1,length(coord3));  %phi4 = NaN(1,length(coord3));       phi_mu4 = NaN(1,length(coord3));

r_sim1 = NaN(1,length(coord3)); r_sim1_ = NaN(1,length(coord3)); rdot_sim1 = NaN(1,length(coord3)); rdot_sim1_ = NaN(1,length(coord3));  diff1 = NaN(1,length(coord3));
r_sim2 = NaN(1,length(coord3)); r_sim2_ = NaN(1,length(coord3)); rdot_sim2 = NaN(1,length(coord3)); rdot_sim2_ = NaN(1,length(coord3));  diff2 = NaN(1,length(coord3));
r_sim3 = NaN(1,length(coord3)); r_sim3_ = NaN(1,length(coord3)); rdot_sim3 = NaN(1,length(coord3)); rdot_sim3_ = NaN(1,length(coord3));  diff3 = NaN(1,length(coord3));
r_sim4 = NaN(1,length(coord3)); r_sim4_ = NaN(1,length(coord3)); rdot_sim4 = NaN(1,length(coord3)); rdot_sim4_ = NaN(1,length(coord3)); diff4 = NaN(1,length(coord3));

%multi = 8*random('Rayleigh',8,[20,1]); gamma = 0.02*rand(20,1) + 0.02*rand(20,1)*i;

multi = randi(10, [20,1]);

for k = 1:1:length(time)-1  
[H1(k+1),H1_(k+1),r_sim1(k+1),r_sim1_(k+1), rdot_sim1(k+1)] = noisysim2(x1,f1,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,time(k+1)-time(k),time(k),multi);
[H2(k+1),H2_(k+1),r_sim2(k+1),r_sim2_(k+1), rdot_sim2(k+1)] = noisysim2(x2,f2,Gt,M,X,PT,2*GT,2*GR,R,sigma,2,k,z,z_prev,time(k+1)-time(k),time(k),multi);
[H3(k+1),H3_(k+1),r_sim3(k+1),r_sim3_(k+1), rdot_sim3(k+1)] = noisysim2(x3,f3,Gt,M,X,PT,2*GT,2*GR,R,sigma,3,k,z,z_prev,time(k+1)-time(k),time(k),multi);   
[H4(k+1),H4_(k+1),r_sim4(k+1),r_sim4_(k+1), rdot_sim4(k+1)] = noisysim2(x4,f4,Gt,M,X,PT,GT,GR,R,sigma,4,k,z,z_prev,time(k+1)-time(k),time(k),multi);  
end
%[v,         phi_mod, r,          r2,            rdot,         rdot2,           diff]

%[H1(k+1),phi_mu1(k+1),r_sim1(k+1),r_sim1_(k+1),rdot_sim1(k+1),rdot_sim1_(k+1),diff1(k+1)]
%[H2(k+1),phi_mu2(k+1),r_sim2(k+1),r_sim2_(k+1),rdot_sim2(k+1),rdot_sim2_(k+1),diff2(k+1)]
%[H3(k+1),phi_mu3(k+1),r_sim3(k+1),r_sim3_(k+1),rdot_sim3(k+1),rdot_sim3_(k+1),diff3(k+1)]
%[H4(k+1),phi_mu4(k+1),r_sim4(k+1),r_sim4_(k+1),rdot_sim4(k+1),rdot_sim4_(k+1),diff4(k+1)]

% rdot_sim1_ = rdot_sim1 + 0.005*rand(1, length(rdot_sim1));% - 0.0025;
% rdot_sim2_ = rdot_sim2 + 0.005*rand(1, length(rdot_sim2));% - 0.0025;
% rdot_sim3_ = rdot_sim3 + 0.005*rand(1, length(rdot_sim3));% - 0.0025;
% rdot_sim4_ = rdot_sim4 + 0.005*rand(1, length(rdot_sim4));% - 0.0025;

% figure
% subplot(4,1,1),plot(time, r_sim1,'LineWidth',3);title('Simulated Radial Distance r1 from Reader $\#1$ in 3D','interpreter','latex');ylabel('r_1 [m]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,2),plot(time, r_sim2,'LineWidth',3);title('Simulated Radial Distance r2 from Reader $\#2$ in 2D','interpreter','latex');ylabel('r_2 [m]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,3),plot(time, r_sim3,'LineWidth',3);title('Simulated Radial Distance r3 from Reader $\#3$ in 3D','interpreter','latex');ylabel('r_3 [m]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,4),plot(time, r_sim4,'LineWidth',3);title('Simulated Radial Distance r4 from Reader $\#4$ in 3D','interpreter','latex');ylabel('r_4 [m]');xlabel('t [s]');grid on; grid minor;

%
figure
subplot(8,1,1),plot(time, H1,'LineWidth',2);title('Simulated H from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(8,1,2),plot(time, mag2db(H1),'LineWidth',2);title('Simulated H from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(8,1,3),plot(time, H2,'LineWidth',2);title('Simulated H from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(8,1,4),plot(time, mag2db(H2),'LineWidth',2);title('Simulated H from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

subplot(8,1,5),plot(time, H3,'LineWidth',2);title('Simulated H from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(8,1,6),plot(time, mag2db(H3),'LineWidth',2);title('Simulated H from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

subplot(8,1,7),plot(time, H4,'LineWidth',2);title('Simulated H from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(8,1,8),plot(time, mag2db(H4),'LineWidth',2);title('Simulated H from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;

figure
subplot(4,1,1),plot(time, rdot_sim1,'LineWidth',2);title('Simulated Phase from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time, rdot_sim2,'LineWidth',2);title('Simulated Phase from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time, rdot_sim3,'LineWidth',2);title('Simulated Phase from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time, rdot_sim4,'LineWidth',2);title('Simulated Phase from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
%%
figure
subplot(4,1,1),plot(time, r_sim1_,'LineWidth',2);title('Simulated Magnitude from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,2),plot(time, r_sim2_,'LineWidth',2);title('Simulated Magnitude from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,3),plot(time, r_sim3_,'LineWidth',2);title('Simulated Magnitude from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
subplot(4,1,4),plot(time, r_sim4_,'LineWidth',2);title('Simulated Magnitude from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
%%

% 
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

%
time4 = 0:length(magD4)/(length(magD4)*6000):length(magD4)/6000 - length(magD4)/(length(magD4)*6000);
time4 = time4*1.238;

mag1   = magD1(145388:315021);     mag2 = magD2(145388:315021);     mag3 = magD3(295791:640912);     mag4 = magD4(145388:315021);
phase1 = phaseD1(145388:315021); phase2 = phaseD2(145388:315021); phase3 = phaseD3(295791:640912); phase4 = phaseD4(145388:315021);
t1 = time1(145388:315021);           t2 = time2(145388:315021);       t3 = time3(295791:640912);      t4 = time4(145388:315021); 


%save('RF.mat','t1','mag1','phase1','t2','mag2','phase2','t3','mag3','phase3','t4','mag4','phase4');

%[time2_4, magD3_4, magD4_4, magD5_4, magD6_4, magD7_4, phaseD2_4] = filterAvg(magD4, phaseD4);
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
%%
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
%%
figure
subplot(411), plot(time, magD12);
subplot(412), plot(time, magD22);
subplot(413), plot(time, magD32);
subplot(414), plot(time, magD42);
%%
%
figure 
subplot(411), plot(time, magD12, 'LineWidth', 2), title('Reader 1 Magnitude'),xlim([35, 155]);
hold on
plot(time, H1, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
hold on
%plot(time2_1, magD3_1, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','SouthEast')
grid on 
grid minor

subplot(412), plot(time, magD22, 'LineWidth', 2), title('Reader 2 Magnitude'),xlim([35, 155]);
hold on
%plot(time, 2.5*H2-0.00029, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
plot(time, H2, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
hold on
%plot(time2_2, magD4_2, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthEast')
grid on 
grid minor

subplot(413), plot(time, magD32, 'LineWidth', 2), title('Reader 3 Magnitude'),xlim([35, 155]);
hold on
%plot(time, 8*H3-0.0011, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
plot(time, H3, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');

hold on
%plot(time2_3, magD4_3, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthEast')
grid on 
grid minor

subplot(414), plot(time, magD42, 'LineWidth', 2), title('Reader 4 Magnitude'),xlim([35, 155]);
hold on
plot(time, H4, 'LineWidth', 2),xlim([35, 155]);
hold on
%plot(time2_4, magD4_4, 'LineWidth', 2),
legend('Measurement','Simulation','Filtered','location','NorthEast'), ylabel('Magnitude [V]');
grid on 
grid minor
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


%%
figure
subplot(811), plot(time1+13.3287, phaseD1,'LineWidth',2), xlim([14, 145])
subplot(812), plot(time, rdot_sim1,'LineWidth',2), xlim([14, 145])
subplot(813), plot(time2+13.3287, phaseD2,'LineWidth',2), xlim([14, 145])
subplot(814), plot(time, rdot_sim2,'LineWidth',2), xlim([14, 145])

subplot(815), plot(time3+13.3287, phaseD3,'LineWidth',2), xlim([14, 145])
subplot(816), plot(time, rdot_sim3,'LineWidth',2), xlim([14, 145])
subplot(817), plot(time4+13.3287, phaseD4,'LineWidth',2), xlim([14, 145])
subplot(818), plot(time, rdot_sim4,'LineWidth',2), xlim([14, 145])
%%
% Sim 2

%
% H12 = NaN(1,length(coord3));     %phi1 = NaN(1,length(coord3));       phi_mu1 = NaN(1,length(coord3));
% H22 = NaN(1,length(coord3));     %phi2 = NaN(1,length(coord3));       phi_mu2 = NaN(1,length(coord3));
% H32 = NaN(1,length(coord3));     %phi3 = NaN(1,length(coord3));       phi_mu3 = NaN(1,length(coord3));
% H42 = NaN(1,length(coord3));     %phi4 = NaN(1,length(coord3));       phi_mu4 = NaN(1,length(coord3));
% 
% r_sim12 = NaN(1,length(coord3)); r_sim12_ = NaN(1,length(coord3)); %rdot_sim1 = NaN(1,length(coord3)); rdot_sim1_ = NaN(1,length(coord3));  diff1 = NaN(1,length(coord3));
% r_sim22 = NaN(1,length(coord3)); r_sim22_ = NaN(1,length(coord3)); %rdot_sim2 = NaN(1,length(coord3)); rdot_sim2_ = NaN(1,length(coord3));  diff2 = NaN(1,length(coord3));
% r_sim32 = NaN(1,length(coord3)); r_sim32_ = NaN(1,length(coord3)); %rdot_sim3 = NaN(1,length(coord3)); rdot_sim3_ = NaN(1,length(coord3));  diff3 = NaN(1,length(coord3));
% r_sim42 = NaN(1,length(coord3)); r_sim42_ = NaN(1,length(coord3)); %rdot_sim4 = NaN(1,length(coord3)); rdot_sim4_ = NaN(1,length(coord3)); diff4 = NaN(1,length(coord3));
% 
% for k = 1:1:length(time)-1  
% [H12(k+1),r_sim12(k+1),r_sim12_(k+1)] = noisysim2(x1,f1,Gt,M,X,PT,GT,GR,R,sigma,1,k,z,z_prev,time(k+1)-time(k));
% [H22(k+1),r_sim22(k+1),r_sim22_(k+1)] = noisysim2(x2,f2,Gt,M,X,PT,2*GT,2*GR,R,sigma,2,k,z,z_prev,time(k+1)-time(k));
% [H32(k+1),r_sim32(k+1),r_sim32_(k+1)] = noisysim(x3,f3,Gt,M,X,PT,2*GT,2*GR,R,sigma,3,k,z,z_prev,time(k+1)-time(k));   
% [H42(k+1),r_sim42(k+1),r_sim42_(k+1)] = noisysim2(x4,f4,Gt,M,X,PT,GT,GR,R,sigma,4,k,z,z_prev,time(k+1)-time(k));  
% end
% 
% figure
% subplot(4,1,1),plot(time, H12,'LineWidth',2);title('Simulated H from Reader $\#1$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,2),plot(time, H22,'LineWidth',2);title('Simulated H from Reader $\#2$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,3),plot(time, H32,'LineWidth',2);title('Simulated H from Reader $\#3$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% subplot(4,1,4),plot(time, H42,'LineWidth',2);title('Simulated H from Reader $\#4$ in 3D','interpreter','latex');ylabel('Magnitude [V]');xlabel('t [s]');grid on; grid minor;
% 
% 
% %%
% figure 
% subplot(411), plot(time1+13.3287, magD1, 'LineWidth', 2), title('Reader 1 Magnitude'),xlim([35, 155]);
% hold on
% plot(time, H12, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
% hold on
% %plot(time2_1, magD3_1, 'LineWidth', 2),
% legend('Measurement','Simulation','Filtered','location','SouthEast')
% grid on 
% grid minor
% 
% subplot(412), plot(time2+13.3287, magD2, 'LineWidth', 2), title('Reader 2 Magnitude'),xlim([35, 155]);
% hold on
% %plot(time, 2.5*H2-0.00029, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
% plot(time, H22, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
% hold on
% %plot(time2_2, magD4_2, 'LineWidth', 2),
% legend('Measurement','Simulation','Filtered','location','NorthEast')
% grid on 
% grid minor
% 
% subplot(413), plot(time3+13.3287, magD3, 'LineWidth', 2), title('Reader 3 Magnitude'),xlim([35, 155]);
% hold on
% %plot(time, 8*H3-0.0011, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
% plot(time, H32, 'LineWidth', 2),xlim([35, 155]);ylabel('Magnitude [V]');
% 
% hold on
% %plot(time2_3, magD4_3, 'LineWidth', 2),
% legend('Measurement','Simulation','Filtered','location','NorthEast')
% grid on 
% grid minor
% 
% subplot(414), plot(time4+13.3287, magD4, 'LineWidth', 2), title('Reader 4 Magnitude'),xlim([35, 155]);
% hold on
% plot(time, H42-0.00001, 'LineWidth', 2),xlim([35, 155]);
% hold on
% %plot(time2_4, magD4_4, 'LineWidth', 2),
% legend('Measurement','Simulation','Filtered','location','NorthEast'), ylabel('Magnitude [V]');
% grid on 
% grid minor