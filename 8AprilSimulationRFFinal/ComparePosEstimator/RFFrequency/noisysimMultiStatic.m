function [v, v2, rmag, rgt, rmeas, rphase, rdot, rdotgt, phi_mod_gt, phi_mod_gt2, phi_mod, phi_mod2, phi_mod3, phi_mod4, offset1, offset2, offset3, offset4, offset5, offset6] = noisysimMultiStatic(x,f,Gt,M,X,PT,GT,GR,R,sigma,k,z,z_prev,phi_prev_mod_gt_load, phi_prev_mod_load, phi_prev_mod4_load, T, Hmeas, offset1, offset2, offset3, offset4, offset5, offset6)

% ============================================= Magnitude =================================================
Delta_f = 0.1*10^9;    lambda = 3*10^8/f;     lambda2 = 3*10^8/(f+Delta_f); lambda3 = 3*10^8/(f-Delta_f); lambda4 = 3*10^8/(f-2*Delta_f);

xcoord = z(1,k)-x(1);  ycoord = z(2,k)-x(2);   zcoord = z(3,k)-x(3);
% ============================================= Multi-path =================================================
% --------------- The dimension of the room is 7.32m * 6.8m * 3m -------------

direct        = sqrt((xcoord)^2+(ycoord)^2+(zcoord)^2);
mul           = 1/direct*exp(-i*2*pi*direct/lambda);
mul2          = 1/direct*exp(-i*2*pi*direct/lambda2);
mul3          = 1/direct*exp(-i*2*pi*direct/lambda3);
mul4          = 1/direct*exp(-i*2*pi*direct/lambda4);

direct_prev   = sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2+(z_prev(3,k)-x(3))^2);
mul_prev      = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda);
mul_prev2     = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda2);
mul_prev3     = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda3);
mul_prev4     = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda4);

drefl1        = sqrt((x(1) + 3)^2   + ((x(2) - z(2,k))/2)^2  + ((x(3) - z(3,k))/2)^2)           + sqrt((z(1,k) + 3)^2 + ((x(2) - z(2,k))/2)^2 + ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef1       = 1/drefl1*exp(-i*2*pi*drefl1/lambda);
mulRef12      = 1/drefl1*exp(-i*2*pi*drefl1/lambda2);
mulRef13      = 1/drefl1*exp(-i*2*pi*drefl1/lambda3);
mulRef14      = 1/drefl1*exp(-i*2*pi*drefl1/lambda4);

drefl1_prev   = sqrt((x(1) + 3)^2   + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2)   + sqrt((z_prev(1,k) + 3)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef1_prev  = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda);
mulRef12_prev = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda2);
mulRef13_prev = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda3);
mulRef14_prev = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda4);

drefl2        = sqrt((4.32- x(1))^2   + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2)     + sqrt((4.32 - z(1,k))^2 + ((x(2) - z(2,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef2       = 1/drefl2*exp(-i*2*pi*drefl2/lambda);
mulRef22      = 1/drefl2*exp(-i*2*pi*drefl2/lambda2);
mulRef23      = 1/drefl2*exp(-i*2*pi*drefl2/lambda3);
mulRef24      = 1/drefl2*exp(-i*2*pi*drefl2/lambda4);

drefl2_prev   = sqrt((4.32- x(1))^2   + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((4.32 - z_prev(1,k))^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef2_prev  = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda);
mulRef22_prev = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda2);
mulRef23_prev = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda3);
mulRef24_prev = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda4);

drefl3        = sqrt((x(2) + 3.36)^2   + ((x(1) - z(1,k))/2)^2     + ((x(3) - z(3,k))/2)^2)      + sqrt((z(2,k) + 3.36)^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef3       = 1/drefl3*exp(-i*2*pi*drefl3/lambda);
mulRef32      = 1/drefl3*exp(-i*2*pi*drefl3/lambda2);
mulRef33      = 1/drefl3*exp(-i*2*pi*drefl3/lambda3);
mulRef34      = 1/drefl3*exp(-i*2*pi*drefl3/lambda4);

drefl3_prev   = sqrt((x(2) + 3.36)^2   + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((z_prev(2,k) + 3.36)^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef3_prev  = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda);
mulRef32_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda2);
mulRef33_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda3);
mulRef34_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda4);

drefl4        = sqrt((3.44 - x(2))^2   + ((x(1) - z(1,k))/2)^2     + ((x(3) - z(3,k))/2)^2)      + sqrt((3.44 - z(2,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef4       = 1/drefl4*exp(-i*2*pi*drefl4/lambda);
mulRef42      = 1/drefl4*exp(-i*2*pi*drefl4/lambda2);
mulRef43      = 1/drefl4*exp(-i*2*pi*drefl4/lambda3);
mulRef44      = 1/drefl4*exp(-i*2*pi*drefl4/lambda4);

drefl4_prev   = sqrt((3.44 - x(2))^2   + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((3.44 - z_prev(2,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef4_prev  = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda);
mulRef42_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda2);
mulRef43_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda3);
mulRef44_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda4);

drefl5        = sqrt((x(3) + 1)^2   + ((x(2) - z(2,k))/2)^2     + ((x(1) - z(1,k))/2)^2)      + sqrt((z(3,k) + 1)^2 + ((x(2) - z(2,k))/2)^2+ ((x(1) - z(1,k))/2)^2);%direct + lambda*50;
mulRef5       = 1/drefl5*exp(-i*2*pi*drefl5/lambda);
mulRef52      = 1/drefl5*exp(-i*2*pi*drefl5/lambda2);
mulRef53      = 1/drefl5*exp(-i*2*pi*drefl5/lambda3);
mulRef54      = 1/drefl5*exp(-i*2*pi*drefl5/lambda4);

drefl5_prev   = sqrt((x(3) + 1)^2   + ((x(2) - z_prev(2,k))/2)^2+ ((x(1) - z_prev(1,k))/2)^2) + sqrt((z_prev(3,k) + 1)^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(1) - z_prev(1,k))/2)^2);%direct + lambda*50;
mulRef5_prev  = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda);
mulRef52_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda2);
mulRef53_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda3);
mulRef54_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda4);

drefl6        = sqrt((2 - x(3))^2  + ((x(1) - z(1,k))/2)^2     + ((x(2) - z(2,k))/2)^2)       + sqrt((2 - z(3,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(2) - z(2,k))/2)^2);%direct + lambda*50;
mulRef6       = 1/drefl6*exp(-i*2*pi*drefl6/lambda);
mulRef62      = 1/drefl6*exp(-i*2*pi*drefl6/lambda2);
mulRef63      = 1/drefl6*exp(-i*2*pi*drefl6/lambda3);
mulRef64      = 1/drefl6*exp(-i*2*pi*drefl6/lambda4);

drefl6_prev   = sqrt((2 - x(3))^2  + ((x(1) - z_prev(1,k))/2)^2 + ((x(2) - z_prev(2,k))/2)^2) + sqrt((2 - z_prev(3,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(2) - z_prev(2,k))/2)^2);%direct + lambda*50;
mulRef6_prev  = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda);
mulRef62_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda2);
mulRef63_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda3);
mulRef64_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda4);

H  = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul+ 0.1*mulRef1 + 0.1*mulRef2 + 0.1*mulRef3 + 0.1*mulRef4 + 0.1*mulRef5 +0.1*mulRef6))^4)); %  + 0.6*mulSum3 %mul+ 0.25* mulSum + 0.25*mulSum2
H2 = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul))^4));

% ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++

H = H^2/abs(H);
%pd = makedist('Rician','s',sqrt(H),'sigma',0.10*sqrt(H));
%H = random(pd).^2;

v = H;
v2 = H2;

rmag  = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*H^2))^(-1/4);
rgt   = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*H2^2))^(-1/4);
rmeas = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*Hmeas^2))^(-1/4);

% ============================================== Phase ====================================================

phi_prev_mod_gt  = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod_gt       = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

phi_prev_mod2_gt = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul_prev2)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod_gt2      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul2)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

phi_prev_mod = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul_prev + 0.1*mulRef1_prev+ 0.1*mulRef2_prev   + 0.1*mulRef3_prev + 0.1*mulRef4_prev + 0.1*mulRef5_prev + 0.1*mulRef6_prev)^4));
phi_mod      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul      + 0.1*mulRef1     + 0.1*mulRef2        + 0.1*mulRef3      + 0.1*mulRef4      + 0.1*mulRef5      + 0.1*mulRef6)^4));

phi_prev_mod2= angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul_prev2+ 0.1*mulRef12_prev+ 0.1*mulRef22_prev+ 0.1*mulRef32_prev+ 0.1*mulRef42_prev+ 0.1*mulRef52_prev+ 0.1*mulRef62_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod2     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul2     + 0.1*mulRef12     + 0.1*mulRef22     + 0.1*mulRef32     + 0.1*mulRef42     + 0.1*mulRef52     + 0.1*mulRef62)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

phi_prev_mod3= angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda3^4*X^2*M)/(4*pi)^4*(mul_prev3+ 0.1*mulRef13_prev+ 0.1*mulRef23_prev+ 0.1*mulRef33_prev+ 0.1*mulRef43_prev+ 0.1*mulRef53_prev+ 0.1*mulRef63_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod3     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda3^4*X^2*M)/(4*pi)^4*(mul3     + 0.1*mulRef13     + 0.1*mulRef23     + 0.1*mulRef33     + 0.1*mulRef43     + 0.1*mulRef53     + 0.1*mulRef63)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

phi_prev_mod4= angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda4^4*X^2*M)/(4*pi)^4*(mul_prev4+ 0.1*mulRef14_prev+ 0.1*mulRef24_prev+ 0.1*mulRef34_prev+ 0.1*mulRef44_prev+ 0.1*mulRef54_prev+ 0.1*mulRef64_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod4     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda4^4*X^2*M)/(4*pi)^4*(mul4     + 0.1*mulRef14     + 0.1*mulRef24     + 0.1*mulRef34     + 0.1*mulRef44     + 0.1*mulRef54     + 0.1*mulRef64)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

if abs(phi_prev_mod_gt - phi_mod_gt) > 2.8
    %offset1 = offset1 + pi*sign(phi_prev_mod_gt - phi_mod_gt);
    offset1 = offset1 + phi_prev_mod_gt - phi_mod_gt;
end

if abs(phi_prev_mod2_gt - phi_mod_gt2) > 2.8
    %offset2 = offset2 + pi*sign(phi_prev_mod2_gt - phi_mod_gt2);
    offset2 = offset2 + phi_prev_mod2_gt - phi_mod_gt2;
end

if abs(phi_prev_mod - phi_mod) > 2.8
    %offset3 = offset3 + pi*sign(phi_prev_mod - phi_mod);
    offset3 = offset3 + phi_prev_mod - phi_mod;
end

if abs(phi_prev_mod2 - phi_mod2) > 2.8
    %offset4 = offset4 + pi*sign(phi_prev_mod2 - phi_mod2);
    offset4 = offset4 + phi_prev_mod2 - phi_mod2;
end

if abs(phi_prev_mod3 - phi_mod3) > 2.8
    %offset5 = offset5 + pi*sign(phi_prev_mod3 - phi_mod3);
    offset5 = offset5 + phi_prev_mod3 - phi_mod3;
end

if abs(phi_prev_mod4 - phi_mod4) > 2.8
    %offset6 = offset6 + pi*sign(phi_prev_mod4 - phi_mod4);
    offset6 = offset6 + phi_prev_mod4 - phi_mod4;
end

phi_mod_gt  = phi_mod_gt + offset1;
phi_mod_gt2 = phi_mod_gt2+ offset2;

phi_mod  = phi_mod + offset3;
phi_mod2 = phi_mod2+ offset4;

phi_mod3 = phi_mod3+ offset5;
phi_mod4 = phi_mod4+ offset6;


diff = phi_mod - phi_prev_mod_load;

%-------------------
%diff = phi_mod - phi_prev_mod;

delta_phi = diff;

% +++++++++++++++++++++++++++++++++++++++++++++ Noise of Phase +++++++++++++++++++++++++++++++++++++++++++++++

new_sigma = sigma/phi_mod_gt + sigma/phi_mod_gt;
phi_noise = 0.04*new_sigma.*rand;

new_sigma2 = sigma/phi_mod_gt2 + sigma/phi_mod_gt2;
phi_noise2 = 0.04*new_sigma2.*rand;

% phi_mod_gt  = exp(phi_noise)*phi_mod_gt;
% phi_mod_gt2 = exp(phi_noise2)*phi_mod_gt2;

% rdot2 = lambda/(4*pi)*delta_phi2*1/T;
if isnan(phi_prev_mod_load)
    phi_prev_mod_load = phi_mod2;
end

% if isnan(phi_prev_mod4_load)
%     phi_prev_mod4_load = phi_mod3;
% end

% phi_mod = phi_mod
% phi_mod2 = phi_mod2
% phi_mod4 = phi_mod4
% phi_mod3 = phi_mod3

rphase_prev = 3*10^8/(4*pi)*(phi_prev_mod_gt - phi_prev_mod2_gt)/Delta_f;
xx = 1:1:4;
ans = polyfit(xx, [phi_mod, phi_mod2, phi_mod3, phi_mod4], 1);

%rphase  = 3*10^8/(4*pi)*1/2*((phi_prev_mod_load - phi_mod2)+(phi_prev_mod4_load - phi_mod3))/Delta_f;%lambda/(4*pi)*phi_mod*1/T; %+ (phi_mod2 - phi_mod4)
%rphase  = 3*10^8/(4*pi)*1/1*((phi_mod - phi_mod2)+0*(phi_mod4 - phi_mod3))/Delta_f;%lambda/(4*pi)*phi_mod*1/T; %+ (phi_mod2 - phi_mod4)
rphase  = 3*10^8/(4*pi)*(ans(1))/Delta_f;%lambda/(4*pi)*phi_mod*1/T; %+ (phi_mod2 - phi_mod4)

rdot   = lambda/(4*pi)*(phi_mod_gt - phi_prev_mod_gt_load)*1/T;
rdotgt = lambda/(4*pi)*diff*1/T;

if abs(rdotgt) > 2
    rdotgt = 0;
end
% if abs(rphase) > 7
%     rphase = 0;
% end

% ============================================= Magnitude =================================================
% Delta_f = 0.5*10^5;    lambda = 3*10^8/f;     lambda2 = 3*10^8/(f+Delta_f); lambda3 = 3*10^8/(f-3*Delta_f); lambda4 = 3*10^8/(f-4*Delta_f);
% 
% xcoord = z(1,k)-x(1);  ycoord = z(2,k)-x(2);   zcoord = z(3,k)-x(3);
% % ============================================= Multi-path =================================================
% % --------------- The dimension of the room is 7.32m * 6.8m * 3m -------------
% 
% direct        = sqrt((xcoord)^2+(ycoord)^2+(zcoord)^2);
% mul           = 1/direct*exp(-i*2*pi*direct/lambda);
% mul2          = 1/direct*exp(-i*2*pi*direct/lambda2);
% mul3          = 1/direct*exp(-i*2*pi*direct/lambda3);
% mul4          = 1/direct*exp(-i*2*pi*direct/lambda4);
% 
% direct_prev   = sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2+(z_prev(3,k)-x(3))^2);
% mul_prev      = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda);
% mul_prev2     = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda2);
% mul_prev3     = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda3);
% mul_prev4     = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda4);
% 
% drefl1        = sqrt((x(1) + 3)^2   + ((x(2) - z(2,k))/2)^2  + ((x(3) - z(3,k))/2)^2)           + sqrt((z(1,k) + 3)^2 + ((x(2) - z(2,k))/2)^2 + ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
% mulRef1       = 1/drefl1*exp(-i*2*pi*drefl1/lambda);
% mulRef12      = 1/drefl1*exp(-i*2*pi*drefl1/lambda2);
% mulRef13      = 1/drefl1*exp(-i*2*pi*drefl1/lambda3);
% mulRef14      = 1/drefl1*exp(-i*2*pi*drefl1/lambda4);
% 
% drefl1_prev   = sqrt((x(1) + 3)^2   + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2)   + sqrt((z_prev(1,k) + 3)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
% mulRef1_prev  = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda);
% mulRef12_prev = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda2);
% mulRef13_prev = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda3);
% mulRef14_prev = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda4);
% 
% drefl2        = sqrt((4.32- x(1))^2   + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2)     + sqrt((4.32 - z(1,k))^2 + ((x(2) - z(2,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
% mulRef2       = 1/drefl2*exp(-i*2*pi*drefl2/lambda);
% mulRef22      = 1/drefl2*exp(-i*2*pi*drefl2/lambda2);
% mulRef23      = 1/drefl2*exp(-i*2*pi*drefl2/lambda3);
% mulRef24      = 1/drefl2*exp(-i*2*pi*drefl2/lambda4);
% 
% drefl2_prev   = sqrt((4.32- x(1))^2   + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((4.32 - z_prev(1,k))^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
% mulRef2_prev  = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda);
% mulRef22_prev = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda2);
% mulRef23_prev = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda3);
% mulRef24_prev = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda4);
% 
% drefl3        = sqrt((x(2) + 3.36)^2   + ((x(1) - z(1,k))/2)^2     + ((x(3) - z(3,k))/2)^2)      + sqrt((z(2,k) + 3.36)^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
% mulRef3       = 1/drefl3*exp(-i*2*pi*drefl3/lambda);
% mulRef32      = 1/drefl3*exp(-i*2*pi*drefl3/lambda2);
% mulRef33      = 1/drefl3*exp(-i*2*pi*drefl3/lambda3);
% mulRef34      = 1/drefl3*exp(-i*2*pi*drefl3/lambda4);
% 
% drefl3_prev   = sqrt((x(2) + 3.36)^2   + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((z_prev(2,k) + 3.36)^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
% mulRef3_prev  = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda);
% mulRef32_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda2);
% mulRef33_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda3);
% mulRef34_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda4);
% 
% drefl4        = sqrt((3.44 - x(2))^2   + ((x(1) - z(1,k))/2)^2     + ((x(3) - z(3,k))/2)^2)      + sqrt((3.44 - z(2,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
% mulRef4       = 1/drefl4*exp(-i*2*pi*drefl4/lambda);
% mulRef42      = 1/drefl4*exp(-i*2*pi*drefl4/lambda2);
% mulRef43      = 1/drefl4*exp(-i*2*pi*drefl4/lambda3);
% mulRef44      = 1/drefl4*exp(-i*2*pi*drefl4/lambda4);
% 
% drefl4_prev   = sqrt((3.44 - x(2))^2   + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((3.44 - z_prev(2,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
% mulRef4_prev  = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda);
% mulRef42_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda2);
% mulRef43_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda3);
% mulRef44_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda4);
% 
% drefl5        = sqrt((x(3) + 1)^2   + ((x(2) - z(2,k))/2)^2     + ((x(1) - z(1,k))/2)^2)      + sqrt((z(3,k) + 1)^2 + ((x(2) - z(2,k))/2)^2+ ((x(1) - z(1,k))/2)^2);%direct + lambda*50;
% mulRef5       = 1/drefl5*exp(-i*2*pi*drefl5/lambda);
% mulRef52      = 1/drefl5*exp(-i*2*pi*drefl5/lambda2);
% mulRef53      = 1/drefl5*exp(-i*2*pi*drefl5/lambda3);
% mulRef54      = 1/drefl5*exp(-i*2*pi*drefl5/lambda4);
% 
% drefl5_prev   = sqrt((x(3) + 1)^2   + ((x(2) - z_prev(2,k))/2)^2+ ((x(1) - z_prev(1,k))/2)^2) + sqrt((z_prev(3,k) + 1)^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(1) - z_prev(1,k))/2)^2);%direct + lambda*50;
% mulRef5_prev  = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda);
% mulRef52_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda2);
% mulRef53_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda3);
% mulRef54_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda4);
% 
% drefl6        = sqrt((2 - x(3))^2  + ((x(1) - z(1,k))/2)^2     + ((x(2) - z(2,k))/2)^2)       + sqrt((2 - z(3,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(2) - z(2,k))/2)^2);%direct + lambda*50;
% mulRef6       = 1/drefl6*exp(-i*2*pi*drefl6/lambda);
% mulRef62      = 1/drefl6*exp(-i*2*pi*drefl6/lambda2);
% mulRef63      = 1/drefl6*exp(-i*2*pi*drefl6/lambda3);
% mulRef64      = 1/drefl6*exp(-i*2*pi*drefl6/lambda4);
% 
% drefl6_prev   = sqrt((2 - x(3))^2  + ((x(1) - z_prev(1,k))/2)^2 + ((x(2) - z_prev(2,k))/2)^2) + sqrt((2 - z_prev(3,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(2) - z_prev(2,k))/2)^2);%direct + lambda*50;
% mulRef6_prev  = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda);
% mulRef62_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda2);
% mulRef63_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda3);
% mulRef64_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda4);
% 
% H  = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul+ 0.1*mulRef1 + 0.1*mulRef2 + 0.1*mulRef3 + 0.1*mulRef4 + 0.1*mulRef5 +0.1*mulRef6))^4)); %  + 0.6*mulSum3 %mul+ 0.25* mulSum + 0.25*mulSum2
% H2 = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul))^4));
% 
% % ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++
% 
% H = H^2/abs(H);
% pd = makedist('Rician','s',sqrt(H),'sigma',0.1*sqrt(H));
% H = random(pd).^2;
% 
% v = H;
% v2 = H2;
% 
% rmag  = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*H^2))^(-1/4);
% rgt   = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*H2^2))^(-1/4);
% rmeas = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*Hmeas^2))^(-1/4);
% 
% % ============================================== Phase ====================================================
% 
% phi_prev_mod_gt  = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
% phi_mod_gt       = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% 
% phi_prev_mod2_gt = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul_prev2)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
% phi_mod_gt2      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul2)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% 
% phi_prev_mod = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul_prev + 0.1*mulRef1_prev+ 0.1*mulRef2_prev+ 0.1*mulRef3_prev+ 0.1*mulRef4_prev+ 0.1*mulRef5_prev+ 0.1*mulRef6_prev)^4));
% phi_mod      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul      + 0.1*mulRef1     + 0.1*mulRef2     + 0.1*mulRef3     + 0.1*mulRef4     + 0.1*mulRef5     + 0.1*mulRef6)^4));
% 
% phi_prev_mod2= angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul_prev2+ 0.1*mulRef12_prev+ 0.1*mulRef22_prev+ 0.1*mulRef32_prev+ 0.1*mulRef42_prev+ 0.1*mulRef52_prev+ 0.1*mulRef62_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
% phi_mod2     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul2     + 0.1*mulRef12     + 0.1*mulRef22     + 0.1*mulRef32     + 0.1*mulRef42     + 0.1*mulRef52     + 0.1*mulRef62)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% 
% phi_prev_mod3= angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda3^4*X^2*M)/(4*pi)^4*(mul_prev3+ 0.1*mulRef13_prev+ 0.1*mulRef23_prev+ 0.1*mulRef33_prev+ 0.1*mulRef43_prev+ 0.1*mulRef53_prev+ 0.1*mulRef63_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
% phi_mod3     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda3^4*X^2*M)/(4*pi)^4*(mul3     + 0.1*mulRef13     + 0.1*mulRef23     + 0.1*mulRef33     + 0.1*mulRef43     + 0.1*mulRef53     + 0.1*mulRef63)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% 
% phi_prev_mod4= angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda4^4*X^2*M)/(4*pi)^4*(mul_prev4+ 0.1*mulRef14_prev+ 0.1*mulRef24_prev+ 0.1*mulRef34_prev+ 0.1*mulRef44_prev+ 0.1*mulRef54_prev+ 0.1*mulRef64_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
% phi_mod4     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda4^4*X^2*M)/(4*pi)^4*(mul4     + 0.1*mulRef14     + 0.1*mulRef24     + 0.1*mulRef34     + 0.1*mulRef44     + 0.1*mulRef54     + 0.1*mulRef64)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% 
% % phi_prev_mod5= angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda4^4*X^2*M)/(4*pi)^4*(mul_prev5+ 0.2*mulRef15_prev+ 0.2*mulRef25_prev+ 0.2*mulRef35_prev+ 0.2*mulRef45_prev+ 0.2*mulRef55_prev+ 0.2*mulRef65_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
% % phi_mod5     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda4^4*X^2*M)/(4*pi)^4*(mul5     + 0.2*mulRef15     + 0.2*mulRef25     + 0.2*mulRef35     + 0.2*mulRef45     + 0.2*mulRef55     + 0.2*mulRef65)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% 
% 
% if abs(phi_prev_mod_gt - phi_mod_gt) > 1
%     offset1 = offset1 + pi*sign(phi_prev_mod_gt - phi_mod_gt);
% end
% 
% if abs(phi_prev_mod2_gt - phi_mod_gt2) > 1
%     offset2 = offset2 + pi*sign(phi_prev_mod2_gt - phi_mod_gt2);
% end
% 
% if abs(phi_prev_mod - phi_mod) > 1
%     offset3 = offset3 + pi*sign(phi_prev_mod - phi_mod);
% end
% 
% if abs(phi_prev_mod2 - phi_mod2) > 1
%     offset4 = offset4 + pi*sign(phi_prev_mod2 - phi_mod2);
% end
% 
% if abs(phi_prev_mod3 - phi_mod3) > 1
%     offset5 = offset5 + pi*sign(phi_prev_mod3 - phi_mod3);
% end
% 
% if abs(phi_prev_mod4 - phi_mod4) > 1
%     offset6 = offset6 + pi*sign(phi_prev_mod4 - phi_mod4);
% end
% 
% phi_mod_gt  = phi_mod_gt + offset1;
% phi_mod_gt2 = phi_mod_gt2+ offset2;
% 
% phi_mod  = phi_mod + offset3;
% phi_mod2 = phi_mod2+ offset4;
% 
% phi_mod3 = phi_mod3+ offset5;
% phi_mod4 = phi_mod4+ offset6;
% 
% diff = phi_mod - phi_prev_mod_load;
% 
% %-------------------
% %diff = phi_mod - phi_prev_mod;
% 
% delta_phi = diff;
% 
% % +++++++++++++++++++++++++++++++++++++++++++++ Noise of Phase +++++++++++++++++++++++++++++++++++++++++++++++
% 
% new_sigma = sigma/phi_mod_gt + sigma/phi_mod_gt;
% phi_noise = 0.004*new_sigma.*rand;
% 
% new_sigma2 = sigma/phi_mod_gt2 + sigma/phi_mod_gt2;
% phi_noise2 = 0.004*new_sigma2.*rand;
% 
% phi_mod_gt  = exp(phi_noise)*phi_mod_gt;
% phi_mod_gt2 = exp(phi_noise2)*phi_mod_gt2;
% 
% % rdot2 = lambda/(4*pi)*delta_phi2*1/T;
% if isnan(phi_prev_mod_load)
%     phi_prev_mod_load = phi_mod2;
% end
% 
% rphase_prev = 3*10^8/(4*pi)*(phi_prev_mod_gt - phi_prev_mod2_gt)/Delta_f;
% rphase  = 3*10^8/(4*pi)*1/2*((phi_mod - phi_mod2)+(phi_mod4 - phi_mod3))/Delta_f;%lambda/(4*pi)*phi_mod*1/T; %+ (phi_mod2 - phi_mod4)
% 
% rdot   = lambda/(4*pi)*(phi_mod_gt - phi_prev_mod_gt_load)*1/T;
% rdotgt = lambda/(4*pi)*diff*1/T;
% 
% if abs(rdotgt) > 2
%     rdotgt = 0;
% end
% if abs(rphase) > 10
%     rphase = 0;
% end

end
