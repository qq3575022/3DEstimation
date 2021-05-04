clc, clear, close all
% cm=1;
% in=0.393*cm
% pr=in*60000%*-1

s=serial('/dev/cu.usbserial-1440','BaudRate',19200,'Parity', 'none', 'DataBits',8, 'StopBits', 1, 'terminator', 'CR');
fopen(s);

x = ones(3,3);        
x(1,1) = 0;       x(1,2) = 0;      x(1,3) = 0;
x(2,1) = 120000;       x(2,2) = 92000;      x(2,3) = 14000;
x(3,1) = 0;  x(3,2) = 0;  x(3,3) = 0;

fprintf(s,'SH XYZ')
i = 3
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/2)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/2)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/2)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;
fprintf(s,dis);      %WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);       %WRITE DISTANCE string to Controller through serial port 
fprintf(s,ac);
fprintf(s,dc);
fprintf(s,'BG XYZ'); %BEGIN Motion along y axis only
