clc, clear, close all
% cm=1;
% in=0.393*cm
% pr=in*60000%*-1

s=serial('/dev/cu.usbserial-1440','BaudRate',19200,'Parity', 'none', 'DataBits',8, 'StopBits', 1, 'terminator', 'CR');
fopen(s);

x = ones(7,3);        
x(1,1) = 0;       x(1,2) = 0;      x(1,3) = 0;
x(2,1) = 10000;   x(2,2) = 0;      x(2,3) = 0;
x(3,1) = 0;       x(3,2) = 0;      x(3,3) = 0;

x(4,1) = 120000;  x(4,2) = 0;      x(4,3) = 0;
x(5,1) = 120000;  x(5,2) = 92000; x(5,3) = 0;
x(6,1) = 120000;  x(6,2) = 92000; x(6,3) = 14000;
x(7,1) = 0;       x(7,2) = 0;      x(7,3) = 0;
x(8,1) = 120000;  x(8,2) = 92000;  x(8,3) = 14000;
x(9,1) = 0;       x(9,2) = 0;      x(9,3) = 0;

pause(8);

%% start

fprintf(s,'SH X')
i = 2;           %Initialize counter
n=6; 

dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))                                            %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/1,(x(i,2)-x(i-1,2))/1,(x(i,3)-x(i-1,3))/1)     %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/30,(x(i,2)-x(i-1,2))/30,(x(i,3)-x(i-1,3))/30) %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/30,(x(i,2)-x(i-1,2))/30,(x(i,3)-x(i-1,3))/30) %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);    %WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);     %WRITE DISTANCE string to Controller through serial port 
%fprintf(s,ac);
%fprintf(s,dc);
fprintf(s,'BG X'); %BEGIN Motion along y axis only
pause(5);
i=i+1;

%% back

fprintf(s,'SH X')
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))                                            %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/1,(x(i,2)-x(i-1,2))/1,(x(i,3)-x(i-1,3))/1)     %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/30,(x(i,2)-x(i-1,2))/30,(x(i,3)-x(i-1,3))/30) %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
%dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/30,(x(i,2)-x(i-1,2))/30,(x(i,3)-x(i-1,3))/30) %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);    %WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);     %WRITE DISTANCE string to Controller through serial port 
%fprintf(s,ac);
%fprintf(s,dc);
fprintf(s,'BG X'); %BEGIN Motion along y axis only
pause(5);
i=i+1;


%% X

fprintf(s,'SH X')
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))                                         %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/2,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/2)  %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)  %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/3,(x(i,3)-x(i-1,3))/3)  %SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);    %WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);     %WRITE DISTANCE string to Controller through serial port 
fprintf(s,ac);
fprintf(s,dc);
fprintf(s,'BG X'); %BEGIN Motion along y axis only
pause(10);
i=i+1;

%% Y

fprintf(s,'SH Y')
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/2,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/2)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/2,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/2)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/2,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/2)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);    %WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);     %WRITE DISTANCE string to Controller through serial port 
fprintf(s,ac);
fprintf(s,dc);
fprintf(s,'BG Y'); %BEGIN Motion along y axis only
pause(10);
i=i+1;

%% Z

fprintf(s,'SH Z')
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/2,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/2,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/2,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);    %WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);     %WRITE DISTANCE string to Controller through serial port 
fprintf(s,ac);
fprintf(s,dc);
fprintf(s,'BG Z'); %BEGIN Motion along y axis only
pause(20);
i=i+1;

%% Back XYZ

fprintf(s,'SH XYZ')
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);      %WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);       %WRITE DISTANCE string to Controller through serial port 
fprintf(s,ac);
fprintf(s,dc);
fprintf(s,'BG XYZ'); %BEGIN Motion along y axis only
pause(20);
i=i+1;

%% Together XYZ

fprintf(s,'SH XYZ')
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
ac =sprintf('AC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input
dc =sprintf('DC %i,%i,%i',(x(i,1)-x(i-1,1))/3,(x(i,2)-x(i-1,2))/2,(x(i,3)-x(i-1,3))/1)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);       %WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);        %WRITE DISTANCE string to Controller through serial port 
fprintf(s,ac);
fprintf(s,dc);
fprintf(s,'BG XYZ');  %BEGIN Motion along y axis only

pause(20);
i=i+1;

%% Back2 XYZ

fprintf(s,'SH XYZ')
dis=sprintf('PA %i,%i,%i',x(i,1),x(i,2),x(i,3))%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

sp =sprintf('SP %i,%i,%i',(x(i,1)-x(i-1,1))/10,(x(i,2)-x(i-1,2))/10,(x(i,3)-x(i-1,3))/10)%SET TRAVEL DISTANCE: Print STRING of SET DISTANCE command from user input

%re = sprintf('TP Y') ;

fprintf(s,dis);       %WRITE DISTANCE string to Controller through serial port 
fprintf(s,sp);        %WRITE DISTANCE string to Controller through serial port 
%fprintf(s,ac);
%fprintf(s,dc);
fprintf(s,'BG XYZ');  %BEGIN Motion along y axis only
