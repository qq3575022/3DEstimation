clc, clear, close all
% ============================  Load instantaneous trilateration r rdot a  ==========================
load('yVectorData.mat');global T; T = mean(diff(tf));
%%
% -----------Load Data Infromation------------
% time: tf    
% radical distances: r1f, r2f, r3f
% radical velocity:  r1dot_e, r2dot_e, r3dot_e
% orientation    psif1
% angular velocity  wf
% acceleration ax1 ay1
% ---------------------------------------------
% Ground truth coordinates of the moving tag path

% ---------------------------------------------
% Get acceleration along x^B y^B
axx = ax1+0.6;  axx = 0.8*axx;  ayy = ay1+0.15; ayy = 0.8*ayy;
AXY = NaN(2,length(axx));       AXY(1,:) = axx; AXY(2,:) = ayy;

% ---------------------------------------------
% Get measurement orientation and angular velocity
offset = 0;
for n = 2:1:length(psif1)
    if psif1(n) - psif1(n-1) < -2
        offset = 2*pi;
    end
    psif1(n) = psif1(n) + offset;  
end
phi = NaN(2,length(axx)); phi(1,:) = psif1 - pi; phi(2,:) = 0.85*wf;

% ---------------------------------------------
% Instantaneous Trilateration
[x_meas, A] = rtoxy(r1f, r2f, r3f, r1dot_e, r2dot_e, r3dot_e, AXY, T);

% ---------------------------------------------
% Get ground-truth
[W_A, W_V, W, XX, rr, AXY_gt, Orient] = loadgtruth(tf);

% ---------------------------------------------
% Get instantaneous trilateration error

% errormeasx = x_meas(1,N+1:length(tf))-XX(1,N+1:length(tf)); 
% errormeasy = x_meas(4,N+1:length(tf))-XX(4,N+1:length(tf)); 
% 
% fprintf('------------- measurement error ------------------')
% error_matricsmeas  = [sqrt(mean(errormeasx)^2+mean(errormeasy)^2), sqrt(var(errormeasx)^2+var(errormeasy)^2),sqrt(var(errormeasx) + var(errormeasy))]
% error_matricsmeasx = [mean(errormeasx), var(errormeasx),sqrt(var(errormeasx))]
% error_matricsmeasy = [mean(errormeasy), var(errormeasy),sqrt(var(errormeasy))]

% Position of readers
x1 = [-0.05, 1.5]; x2 = [2, 3.0]; x3 = [2.7, 0.05];

% Set timestep
time_step    = 0.1;   iternum    = 200; 
time_stepPos = 0.1;   iternumPos = 200;  

% ++++++++++++++++++++++   Input: Length of NLS estimation  ++++++++++++++++
meanrr = NaN(6,40);  stdrr = NaN(6,40);  time = linspace(1,40,40);
meanx = NaN(6,40);   stdx = NaN(6,40);
meany = NaN(6,40);   stdy = NaN(6,40);

%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ==========================  Estimation R Rdot Acc Orientation ======================
N = 5;
y = NaN(N*9,length(tf)-N);

% Get yN
for l = 1:1:length(tf)-N
y(:,l) = getyNPVAO(r1f,r1dot_e,r2f,r2dot_e,r3f,r3dot_e,phi(1,:),AXY(1,:),AXY(2,:),l,N);
end

% Initiate x, e and temproary variable
x = [1.5;0;0;1;0;0;3;0;0]*ones(1,length(tf)-N);e = NaN(1,length(tf)-N);

for m = 1:1:length(tf)-N %time 1 t
    
x(:,m)=lsqnonlin(@(xx)getEPVAO(y(:,m), xx, N),[1.5;0;0;1;0;0;3;0;0]);
end

% Get F^(N-1)*x
x = [1,T,0,0,0,0,0,0,0;
     0,1,T,0,0,0,0,0,0;
     0,0,1,0,0,0,0,0,0;
     0,0,0,1,T,0,0,0,0;
     0,0,0,0,1,T,0,0,0;
     0,0,0,0,0,1,0,0,0;
     0,0,0,0,0,0,1,T,0;
     0,0,0,0,0,0,0,1,T;
     0,0,0,0,0,0,0,0,1]^(N-1)*x;

% Get estimation error
errorestx = x(1,:)-XX(1,N+1:length(tf)); erroresty = x(4,:)-XX(4,N+1:length(tf)); 


%%
N = 5; i = 1; j = 1; k = 1; l = 1;
% Initiate x, e and temproary variable
x = [0.5;0.5;0.5; 0.5;0.5;0.5; 1.1;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]*ones(1,length(time)-N);e = NaN(1,length(time)-N);
factor = 3;

for m = 1:1:length(time)-N
    
    [y, i, j, l, index, len] = gety(r1, rdot1, r2, rdot2, r3+0.05, rdot3, r4, rdot4, angle(1,:), gyro(1,:), angle(2,:), gyro(2,:), angle(3,:), gyro(3,:), acc(1,:), acc(2,:), acc(3,:), accT, gyroT, rtime, time, i, j, l, m, N, factor);

    x(:,m) = lsqnonlin(@(xx)getNLE(y, xx, N, index, len, time, m, factor),[0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5; 0.5;0.5;0.5]);

    F = getF(time, m, m+N);
    x(:,m) = F*x(:,m);

end





%%
figure
subplot(211), plot(tf, wf)
subplot(212), plot(tf, psif1)