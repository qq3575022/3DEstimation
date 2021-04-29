clc, clear, close all
%%
%data=readtable('2.csv','Delimiter', ',');  g=9.7953;
load('RF.mat');
phase11 = phase1(21511:38779); 
time1 = t1(21511:38779);

phaseUwrap1 = NaN(1, length(phase11)); %phaseUwrap2 = NaN(1, length(phase2)); phaseUwrap3 = NaN(1, length(phase3)); phaseUwrap4 = NaN(1, length(phase4));
offset1 = zeros(1, length(phase11)); %offset2 = 0; offset3 = 0; offset4 = 0; 

Tstart = [34.7909, 35.08, 35.3226, 35.4522, 35.6177, 35.7629, 35.8661, 35.9713, 36.0867, 36.202, 36.3986, 36.5241, 36.6302, 36.7583, 36.9314, 37.0814, 37.2597, 37.6106];
Tend = [34.8082, 35.0886, 35.3282, 35.4538, 35.6204, 35.7662, 35.8698, 35.9717, 36.0885, 36.2032, 36.3995, 36.5276, 36.6316, 36.763, 36.9374, 37.0857, 37.2708, 37.6195];

index = 1;
offset = 0;

% for i = 1:1:length(phase11)-1
%     
%     if (abs(Tstart(index) - time1(i))<0.0001)
%         offset1(i) = 
%     end
%     if abs(phase11(i+1) - phase11(i)) > 1
%         offset1 = offset1 - (phase11(i+1) - phase11(i));
%     end
%     
%     phaseUwrap1(i+1) = phase1(i+1) + offset1;
%    
% end
phaseDiff = diff(phase11);

phaseDiff = medfilt1(phaseDiff,300)*1/(2.0634e-03);
% accy = medfilt1(acc_data_y,50) - g;
% accz = medfilt1(acc_data_z,50);

figure
plot(phaseDiff)

for i = 1:1:length(phaseDiff)
    if (abs(phaseDiff(i)) > 1)
        phaseDiff(i) = 0;
    end
end

figure
plot(time1(2:end), phaseDiff*1/(2.0634e-03));
%%
figure
subplot(311), plot(t1(21511:38779), phase11);hold on; plot(t1(21512:38779), diff(phase11));
subplot(312),  plot(t1(21511:38779), phaseUwrap1);
subplot(313), plot(t1(21511:38779), unwrap(phaseUwrap1));


phase33 = phase3(21511:38779);

figure
plot(phase33)
%%
phaseDiff = diff(phase33);

phaseDiff = medfilt1(phaseDiff,300)*1/(2.0634e-03);
% accy = medfilt1(acc_data_y,50) - g;
% accz = medfilt1(acc_data_z,50);

figure
plot(phaseDiff)

for i = 1:1:length(phaseDiff)
    if (abs(phaseDiff(i)) > 1)
        phaseDiff(i) = 0;
    end
end

figure
plot(time1(2:end), phaseDiff*1/(2.0634e-03));
%%
figure
subplot(311), plot(t1(21511:38779), phase11);hold on; plot(t1(21512:38779), diff(phase11));
subplot(312),  plot(t1(21511:38779), phaseUwrap1);
subplot(313), plot(t1(21511:38779), unwrap(phaseUwrap1));

% for i = 1:1:length(phase1)-1
%     if abs(phase1(i+1) - phase1(i)) > 2.5
%         offset1 = offset1 - pi;
%     end
%     
%     if abs(phase2(i+1) - phase2(i)) > 2.5
%         offset2 = offset2 - pi;
%     end
%     
%     if abs(phase3(i+1) - phase3(i)) > 2.5
%         offset3 = offset3 + pi;
%     end
%     
%     if abs(phase4(i+1) - phase4(i)) > 2.5
%         offset4 = offset4 + pi;
%     end
%     
%     phaseUwrap1(i+1) = phase1(i+1) + offset1;
%     phaseUwrap2(i+1) = phase2(i+1) + offset2;
%     phaseUwrap3(i+1) = phase3(i+1) + offset3;
%     phaseUwrap4(i+1) = phase4(i+1) + offset4;
% end