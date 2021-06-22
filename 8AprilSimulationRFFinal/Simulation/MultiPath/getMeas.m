function [magD12, magD22, magD32, magD42, phaseD12, phaseD22, phaseD32, phaseD42] = getMeas(time)

% +++++++++++++++++++++++++++++++++++++++++++ Reader 1 ++++++++++++++++++++++++++++++++++++++
D0_1 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader1_6.bin', 1000000000, 1);
D1 = D0_1;%(85000:445000,1);

magD1 = abs(D1);
phaseD1 = angle(D1);

time1 = 0:length(magD1)/(length(magD1)*6000):length(magD1)/6000 - length(magD1)/(length(magD1)*6000);
time1 = time1*1.238;


% +++++++++++++++++++++++++++++++++++++++++++ Reader 2 ++++++++++++++++++++++++++++++++++++++
D02 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader2_6.bin', 1000000000, 1);
D2 = D02;%(152000:512000,1);

magD2 = abs(D2);
phaseD2 = angle(D2);

time2 = 0:length(magD2)/(length(magD2)*6000):length(magD2)/6000 - length(magD2)/(length(magD2)*6000);
time2 = time2*1.238;


% +++++++++++++++++++++++++++++++++++++++++++ Reader 3 ++++++++++++++++++++++++++++++++++++++
D03 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader3_6.bin', 1000000000, 1);
D3 = D03(31558:end,1);

magD3 = abs(D3);
phaseD3 = angle(D3);


time3 = 0:length(magD3)/(length(magD3)*12000):length(magD3)/12000 - length(magD3)/(length(magD3)*12000);
time3 = time3*1.217;

%[time2_3, magD3_3, magD4_3, magD5_3, magD6_3, magD7_3, phaseD2_3] = filterAvg2(magD3, phaseD3);

% +++++++++++++++++++++++++++++++++++++++++++ Reader 4 ++++++++++++++++++++++++++++++++++++++
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
magD12(1) = magD1(1);         magD22(1) = magD1(1);         magD32(1) = magD1(1);         magD42(1) = magD1(1); 

phaseD12 = NaN(1,length(time)); phaseD22 = NaN(1,length(time)); phaseD32 = NaN(1,length(time)); phaseD42 = NaN(1,length(time));
phaseD12(1) = phase1(1);        phaseD22(1) = phase2(1);        phaseD32(1) = phase3(1);        phaseD42(1) = phase4(1); 

indexT = 2;

for indexMag = 1:1:length(time1)
%     timeT = time(indexT)
%     time1T = time1(indexMag)+13.3287
    if abs(time1(indexMag)+13.3287 - time(indexT)) < 0.0002 || ((time1(indexMag)+13.3287) > time(indexT))
        %indexT
        if indexT  <= length(time)
            magD12(indexT) = magD1(indexMag);
            magD22(indexT) = magD2(indexMag);
            
            phaseD12(indexT) = phaseD1(indexMag);
            phaseD22(indexT) = phaseD2(indexMag);
            
            if indexMag < length(time4)
                magD42(indexT)   = magD4(indexMag);
                phaseD42(indexT) = phaseD4(indexMag);
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
            phaseD32(indexT) = phaseD3(indexMag);
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

end