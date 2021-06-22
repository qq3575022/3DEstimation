clc, clear, close all

% Reader 1
D1 = read_complex_binary ('0304/0304_reader1_single.bin', 1000000000, 1);

magD1 = abs(D1);
phaseD1 = angle(D1);

% 1
D1_2 = read_complex_binary ('0304/0304_reader1.bin', 1000000000, 1);

magD1_2 = abs(D1_2);
phaseD1_2 = angle(D1_2);

% 2
D1_3 = read_complex_binary ('0304/0304_reader1_2.bin', 1000000000, 1);

magD1_3 = abs(D1_3);
phaseD1_3 = angle(D1_3);

% 3
D1_4 = read_complex_binary ('0304/0304_reader1_3.bin', 1000000000, 1);

magD1_4 = abs(D1_4);
phaseD1_4 = angle(D1_4);

figure
subplot(311), plot(magD1), ylim([0, 2*10^(-4)]), grid on, grid minor, title('Single Measurement')
subplot(312), plot(magD1_2), ylim([0, 2*10^(-4)]), grid on, grid minor, title('Combine Measurement 3')
subplot(313), plot(magD1_4), ylim([0, 2*10^(-4)]), grid on, grid minor, title('Combine Measurement 4')

% Reader 2
D2 = read_complex_binary ('0304/0304_reader2_single.bin', 1000000000, 1);

magD2 = abs(D2);
phaseD2 = angle(D2);

% 1
D2_2 = read_complex_binary ('0304/0304_reader2.bin', 1000000000, 1);

magD2_2 = abs(D2_2);
phaseD2_2 = angle(D2_2);

% 2
D2_3 = read_complex_binary ('0304/0304_reader2_2.bin', 1000000000, 1);

magD2_3 = abs(D2_3);
phaseD2_3 = angle(D2_3);

% 3
D2_4 = read_complex_binary ('0304/0304_reader2_3.bin', 1000000000, 1);

magD2_4 = abs(D2_4);
phaseD2_4 = angle(D2_4);

figure
subplot(311), plot(magD2), ylim([0, 3*10^(-4)]), grid on, grid minor, title('Single Measurement')
subplot(312), plot(magD2_2), ylim([0, 3*10^(-4)]), grid on, grid minor, title('Combine Measurement 2')
subplot(313), plot(magD2_4), ylim([0, 2*10^(-4)]), grid on, grid minor, title('Combine Measurement 4')

% Reader 3
D3 = read_complex_binary ('0304/0304_reader3_single.bin', 1000000000, 1);

magD3 = abs(D3);
phaseD3 = angle(D3);

% 1
D3_2 = read_complex_binary ('0304/0304_reader3.bin', 1000000000, 1);

magD3_2 = abs(D3_2);
phaseD3_2 = angle(D3_2);

% 2
D3_3 = read_complex_binary ('0304/0304_reader3_2.bin', 1000000000, 1);

magD3_3 = abs(D3_3);
phaseD3_3 = angle(D3_3);

% 3
D3_4 = read_complex_binary ('0304/0304_reader3_3.bin', 1000000000, 1);

magD3_4 = abs(D3_4);
phaseD3_4 = angle(D3_4);

figure
subplot(311), plot(magD3), ylim([0, 3*10^(-3)]), grid on, grid minor, title('Single Measurement')
subplot(312), plot(magD3_2), ylim([0, 1*10^(-3)]), grid on, grid minor, title('Combine Measurement 2')
subplot(313), plot(magD3_4), ylim([0, 0.7*10^(-3)]), grid on, grid minor, title('Combine Measurement 4')


% Reader 4
D4 = read_complex_binary ('0304/0304_reader4_single.bin', 1000000000, 1);

magD4 = abs(D4);
phaseD4 = angle(D4);


% 1
D4_2 = read_complex_binary ('0304/0304_reader4.bin', 1000000000, 1);

magD4_2 = abs(D4_2);
phaseD4_2 = angle(D4_2);

% 2
D4_3 = read_complex_binary ('0304/0304_reader4_2.bin', 1000000000, 1);

magD4_3 = abs(D4_3);
phaseD4_3 = angle(D4_3);

% 3
D4_4 = read_complex_binary ('0304/0304_reader4_3.bin', 1000000000, 1);

magD4_4 = abs(D4_4);
phaseD4_4 = angle(D4_4);

figure
subplot(311), plot(magD4), ylim([0, 2.5*10^(-4)]), grid on, grid minor, title('Single Measurement')
subplot(312), plot(magD4_2), ylim([0, 2*10^(-4)]), grid on, grid minor, title('Combine Measurement 2')
subplot(313), plot(magD4_4), ylim([0, 1.5*10^(-4)]), grid on, grid minor, title('Combine Measurement 4')



% figure 
% subplot(411), plot(time1, magD1), title('Reader 1')
% hold on
% plot(td3, 0.001*8.22*H1, 'LineWidth', 2), 
% hold on
% plot(time2_1, magD3_1, 'LineWidth', 2),
% legend('Measurement','Simulation','Filtered','location','NorthWest')
% grid on 
% grid minor
% 
% subplot(412), plot(time2, magD2), title('Reader 2')
% hold on
% plot(td3, 0.001*10.61344941*H2, 'LineWidth', 2), 
% hold on
% plot(time2_2, magD4_2, 'LineWidth', 2),
% legend('Measurement','Simulation','Filtered','location','NorthWest')
% grid on 
% grid minor
% 
% subplot(413), plot(time3, magD3), title('Reader 3')
% hold on
% plot(td3, 0.001*42.019839*H3, 'LineWidth', 2), 
% hold on
% plot(time2_3, magD4_3, 'LineWidth', 2),
% legend('Measurement','Simulation','Filtered','location','NorthWest')
% grid on 
% grid minor
% 
% subplot(414), plot(time4, magD4), title('Reader 4')
% hold on
% plot(td3, 0.001*3.779274583*H3, 'LineWidth', 2), 
% hold on
% plot(time2_4, magD4_4, 'LineWidth', 2),
% legend('Measurement','Simulation','Filtered','location','NorthWest')
% grid on 
% grid minor
% 
% % subplot(413), plot(phaseD1), title('Raw Phase Measurement')
% % grid on 
% % grid minor
% % 
% % subplot(414), plot(phaseD2), title('Averaged Phase Measurement')
% % grid on 
% % grid minor
% 
% %
% % mean(d4(1224:1924))