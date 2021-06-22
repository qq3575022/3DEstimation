clc, clear, close all
% Reader 1
%D1 = read_complex_binary ('/Users/Joanna/Documents/MATLAB/3D/MeasurementOct/1005/1005data/1005_reader1_xyz_4.bin', 1000000000, 1);
D1 = read_complex_binary ('0301_test5.bin', 1000000000, 1);

magD1 = abs(D1);
phaseD1 = angle(D1);
%d2=unwrap(phaseD1);

phaseD1T = zeros(1, length(phaseD1));
gtR = zeros(1, length(phaseD1T ));

for k = 1:1:length(phaseD1)
    phaseD1T(k) = k/10000;
    %d3(i-1) = (d2(i) - d2(i-1))/10^(-6);
    if phaseD1T(k) <=9.24 
        gtR(k) = magD1(92400);
    end
    
    if phaseD1T(k) >9.24 && phaseD1T(k) <=12.24
        gtR(k) = magD1(92400) - 0.5*0.04/3*0.001*(phaseD1T(k) - 9.24)^2;
    end
    
    if phaseD1T(k) > 12.24 && phaseD1T(k) <= 19.24

        gtR(k) = gtR(122400) - 0.04*0.001*(phaseD1T(k) - 12.24);
    end
    
    if phaseD1T(k) >=19.24 && phaseD1T(k) <22.24
        gtR(k) = gtR(192400) - 0.04*0.0001*(phaseD1T(k) - 19.24) + 0.5*0.04/3*0.0001*(phaseD1T(k) - 19.24)^2;
    end
    
    if phaseD1T(k) >=22.24 
        gtR(k) = mean(magD1(222400:length(magD1)));
    end
    
end


phaseD2 = zeros(1, round(length(phaseD1)/100));
phaseD2T = zeros(1, length(phaseD2));

for i = 1:1:length(phaseD1)/100
    phaseD2(i) = mean(phaseD1(100*(i-1)+1:100*i));
    phaseD2T(i) = i/100;
    %d3(i-1) = (d2(i) - d2(i-1))/10^(-6);
end

d3 = zeros(1, length(phaseD2) - 1);

% for i = 2:1:length(d3)
%     %phaseD2(i) = mean(phaseD1(1000*(i-1)+1:1000*i));
%     d3(i-1) = (phaseD2(i) - phaseD2(i-1));%/10^(-6);
%     
% end

d3 = diff(phaseD2)/10e-3;

for i = 2:length(d3)
    if (d3(i) - d3(i-1)) > 8
        d3(i) = d3(i-1);
    end
end

d4  = zeros(1, length(d3));
d4T = zeros(1, length(d3));

for j = 1:1:length(d4)-6
    %phaseD2(i) = mean(phaseD1(1000*(i-1)+1:1000*i));
    d4(j) = mean(d3(j:5+j));%/10^(-6);
    d4T(j) = j/100;
    
end

gt = zeros(1, length(d4));

for k = 1:length(phaseD2T)-1
    if phaseD2T(k) >9.24 && phaseD2T(k) <= 12.24
        gt(k) = -5.5643/3*(phaseD2T(k) - 9.24);
    end
    
    if phaseD2T(k) >12.24 && phaseD2T(k) <= 19.24
        gt(k) = -5.5643;
    end
    
    if phaseD2T(k) >19.24 && phaseD2T(k) <=22.24
        gt(k) = -5.5643 + 5.5643/3*(phaseD2T(k) - 19.24);
    end
    
end
% phaseD2 = angle(D2);
% phaseD3 = angle(D3);
%%
figure 
subplot(411), plot(phaseD1T, magD1), title('Raw Magnitude Measurement')
hold on, plot(phaseD1T, gtR, 'LineWidth', 2) 
legend('measurement magnitude', 'groundtruth magnitude')
grid on 
grid minor
subplot(412), plot(phaseD1T, phaseD1), title('Raw Phase Measurement')
grid on 
grid minor
subplot(413), plot(phaseD2T(1:length(phaseD2T)-1), d3),   title('Difference of Phase') 
grid on
grid minor
subplot(414), plot(phaseD2T(1:length(phaseD2T)-1), d4, 'LineWidth', 2), title('Filtered Phase Difference')
hold on, plot(phaseD2T(1:length(phaseD2T)-1), gt, 'LineWidth', 2);
legend('filtered phase', 'velocity')
grid on 
grid minor


%%
mean(d4(1224:1924))