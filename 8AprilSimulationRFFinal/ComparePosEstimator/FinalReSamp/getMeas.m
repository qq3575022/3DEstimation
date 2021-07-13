function [magD1, magD2, mag3, magD4, t1] = getMeas()

% +++++++++++++++++++++++++++++++++++++++++++ Reader 1 ++++++++++++++++++++++++++++++++++++++
D01 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader1_6.bin', 1000000000, 1);
D1 = D01;%(85000:445000,1);

magD1 = abs(D1)+0.0001484;
phaseD1 = angle(D1);

time1 = 0:length(magD1)/(length(magD1)*6000):length(magD1)/6000 - length(magD1)/(length(magD1)*6000);
time1 = time1*1.238;


% +++++++++++++++++++++++++++++++++++++++++++ Reader 2 ++++++++++++++++++++++++++++++++++++++
D02 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader2_6.bin', 1000000000, 1);
D2 = D02;%(152000:512000,1);

magD2 = abs(D2)+0.0005392;
phaseD2 = angle(D2);

time2 = 0:length(magD2)/(length(magD2)*6000):length(magD2)/6000 - length(magD2)/(length(magD2)*6000);
time2 = time2*1.238;


% +++++++++++++++++++++++++++++++++++++++++++ Reader 3 ++++++++++++++++++++++++++++++++++++++
D03 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader3_6.bin', 1000000000, 1);
D3 = D03(31558:end,1);

magD3 = abs(D3)+0.00020778;
phaseD3 = angle(D3);


time3 = 0:length(magD3)/(length(magD3)*12000):length(magD3)/12000 - length(magD3)/(length(magD3)*12000);
time3 = time3*1.217;

% +++++++++++++++++++++++++++++++++++++++++++ Reader 4 ++++++++++++++++++++++++++++++++++++++
D04 = read_complex_binary ('/Users/Joanna/Documents/PhD/Git/Dissertation-Git/3DEstimation/6Data/0422Reader/0422_reader4_6.bin', 1000000000, 1);
D4 = D04;%(150000:510000,1);

magD4 = abs(D4)+0.00015078;
phaseD4 = angle(D4);

time4 = 0:length(magD4)/(length(magD4)*6000):length(magD4)/6000 - length(magD4)/(length(magD4)*6000);
time4 = time4*1.238;

%
mag1   = magD1(145388:315021);     mag2 = magD2(145388:315021);     magD32 = magD3(295791:640912);     mag4 = magD4(145388:315021);
phase1 = phaseD1(145388:315021); phase2 = phaseD2(145388:315021); phase3 = phaseD3(295791:640912); phase4 = phaseD4(145388:315021);
%t1 = time1(145388:315021)'+13.3287;
t2 = time2(145388:315021)+13.3287;t3 = time3(295791:640912)+13.3287;t4 = time4(145388:315021)+13.3287; 

t1 = time1'+13.3287;%t2 = time2+13.3287;t3 = time3+13.3287;t4 = time4+13.3287; 

%
mag3 = NaN(length(t1),1);
mag3(1) = magD1(1);   


indexT = 1;

for indexMag = 1:1:length(time3)
    if indexT  <= length(time1)&& (abs(time3(indexMag) - time1(indexT)) < 0.0002 || ((time3(indexMag)) > time1(indexT)))
        
        mag3(indexT) = magD3(indexMag);
        
        indexT = indexT + 1;
    end
end


end
