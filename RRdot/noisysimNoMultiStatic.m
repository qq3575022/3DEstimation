function [v, v2, rmag, rgt, rmeas, rphase, rdot, rdotgt, phi_mod_gt, phi_mod_gt2, phi_mod, phi_mod2, offset1, offset2, offset3, offset4, offset5, offset6] = noisysimNoMultiStatic(x, f, Gt,X,PT,GT, R,sigma,rice, k,z,z_prev,phi_prev_mod_gt_load, phi_prev_mod_load, T, Hmeas, offset1, offset2, offset3, offset4, offset5, offset6)

% ============================================= Magnitude =================================================
Delta_f = 10^6;    lambda = 3*10^8/f;     lambda2 = 3*10^8/(f+Delta_f); lambda4 = 3*10^8/(f-4*Delta_f); lambda3 = 3*10^8/(f-3*Delta_f); 

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

% Reflection Wall 1
drefl1        = sqrt((x(1) + 3)^2   + ((x(2) - z(2,k))/2)^2  + ((x(3) - z(3,k))/2)^2)           + sqrt((z(1,k) + 3)^2       + ((x(2) - z(2,k))/2)^2 + ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef1       = 1/drefl1*exp(-i*2*pi*drefl1/lambda);
mulRef12      = 1/drefl1*exp(-i*2*pi*drefl1/lambda2);

drefl1_prev   = sqrt((x(1) + 3)^2   + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2)   + sqrt((z_prev(1,k) + 3)^2   + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef1_prev  = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda);
mulRef12_prev = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda2);

% Reflection Wall 2
drefl2        = sqrt((4.32- x(1))^2   + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2)     + sqrt((4.32 - z(1,k))^2     + ((x(2) - z(2,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef2       = 1/drefl2*exp(-i*2*pi*drefl2/lambda);
mulRef22      = 1/drefl2*exp(-i*2*pi*drefl2/lambda2);

drefl2_prev   = sqrt((4.32- x(1))^2   + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((4.32 - z_prev(1,k))^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef2_prev  = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda);
mulRef22_prev = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda2);

% Reflection Wall 3
drefl3        = sqrt((x(2) + 3.36)^2   + ((x(1) - z(1,k))/2)^2     + ((x(3) - z(3,k))/2)^2)      + sqrt((z(2,k) + 3.36)^2      + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef3       = 1/drefl3*exp(-i*2*pi*drefl3/lambda);
mulRef32      = 1/drefl3*exp(-i*2*pi*drefl3/lambda2);

drefl3_prev   = sqrt((x(2) + 3.36)^2   + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((z_prev(2,k) + 3.36)^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef3_prev  = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda);
mulRef32_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda2);

% Reflection Wall 4
drefl4        = sqrt((3.44 - x(2))^2   + ((x(1) - z(1,k))/2)^2     + ((x(3) - z(3,k))/2)^2)      + sqrt((3.44 - z(2,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef4       = 1/drefl4*exp(-i*2*pi*drefl4/lambda);
mulRef42      = 1/drefl4*exp(-i*2*pi*drefl4/lambda2);

drefl4_prev   = sqrt((3.44 - x(2))^2   + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((3.44 - z_prev(2,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef4_prev  = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda);
mulRef42_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda2);

% Reflection Wall 5
drefl5        = sqrt((x(3) + 1)^2   + ((x(2) - z(2,k))/2)^2     + ((x(1) - z(1,k))/2)^2)      + sqrt((z(3,k) + 1)^2 + ((x(2) - z(2,k))/2)^2+ ((x(1) - z(1,k))/2)^2);%direct + lambda*50;
mulRef5       = 1/drefl5*exp(-i*2*pi*drefl5/lambda);
mulRef52      = 1/drefl5*exp(-i*2*pi*drefl5/lambda2);

drefl5_prev   = sqrt((x(3) + 1)^2   + ((x(2) - z_prev(2,k))/2)^2+ ((x(1) - z_prev(1,k))/2)^2) + sqrt((z_prev(3,k) + 1)^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(1) - z_prev(1,k))/2)^2);%direct + lambda*50;
mulRef5_prev  = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda);
mulRef52_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda2);

% Reflection Wall 6
drefl6        = sqrt((2 - x(3))^2  + ((x(1) - z(1,k))/2)^2     + ((x(2) - z(2,k))/2)^2)       + sqrt((2 - z(3,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(2) - z(2,k))/2)^2);%direct + lambda*50;
mulRef6       = 1/drefl6*exp(-i*2*pi*drefl6/lambda);
mulRef62      = 1/drefl6*exp(-i*2*pi*drefl6/lambda2);

drefl6_prev   = sqrt((2 - x(3))^2  + ((x(1) - z_prev(1,k))/2)^2 + ((x(2) - z_prev(2,k))/2)^2) + sqrt((2 - z_prev(3,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(2) - z_prev(2,k))/2)^2);%direct + lambda*50;
mulRef6_prev  = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda);
mulRef62_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda2);

H  = abs(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul+ 0.1*mulRef1 + 0.1*mulRef2 + 0.1*mulRef3 + 0.1*mulRef4 + 0.1*mulRef5 +0.1*mulRef6));
H2 = abs(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*mul);

% ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++

H = H2^2/abs(H2);
pd = makedist('Rician','s',sqrt(H),'sigma',rice*sqrt(H));
H = random(pd).^2;

v = H;
v2 = H2;

rmag  = sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi*H); 
rgt   = sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi*H2); 
rmeas = sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi*Hmeas);

% ============================================== Phase ====================================================

phi_prev_mod_gt  = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul_prev));
phi_mod_gt       = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul));

phi_prev_mod2_gt = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul_prev2));
phi_mod_gt2      = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul2));

phi_prev_mod = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul_prev + 0.1*mulRef1_prev+ 0.1*mulRef2_prev+ 0.1*mulRef3_prev+ 0.1*mulRef4_prev+ 0.1*mulRef5_prev+ 0.1*mulRef6_prev));
phi_mod      = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul      + 0.1*mulRef1     + 0.1*mulRef2     + 0.1*mulRef3     + 0.1*mulRef4     + 0.1*mulRef5     + 0.1*mulRef6));

phi_prev_mod2= angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul_prev2+ 0.1*mulRef12_prev+ 0.1*mulRef22_prev+ 0.1*mulRef32_prev+ 0.1*mulRef42_prev+ 0.1*mulRef52_prev+ 0.1*mulRef62_prev));
phi_mod2     = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul2     + 0.1*mulRef12     + 0.1*mulRef22     + 0.1*mulRef32     + 0.1*mulRef42     + 0.1*mulRef52     + 0.1*mulRef62));

phi_prev_mod3_gt = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul_prev3));
phi_mod_gt3      = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul3));

phi_prev_mod4_gt = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul_prev4));
phi_mod_gt4      = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul4));


if abs(phi_prev_mod_gt - phi_mod_gt) > 3
    offset1 = offset1 + phi_prev_mod_gt - phi_mod_gt;

end

if abs(phi_prev_mod2_gt - phi_mod_gt2) > 3
    offset2 = offset2 + phi_prev_mod2_gt - phi_mod_gt2;
end

if abs(phi_prev_mod - phi_mod) > 3
    offset3 = offset3 + phi_prev_mod - phi_mod;
end

if abs(phi_prev_mod2 - phi_mod2) > 3
    offset4 = offset4 + phi_prev_mod2 - phi_mod2;
end

if abs(phi_prev_mod3_gt - phi_mod_gt3) > 3
    offset5 = offset5 + phi_prev_mod3_gt - phi_mod_gt3;
end

if abs(phi_prev_mod4_gt - phi_mod_gt4) > 3
    offset6 = offset6 + phi_prev_mod4_gt - phi_mod_gt4;
end

phi_mod_gt  = phi_mod_gt + offset1;
phi_mod_gt2 = phi_mod_gt2+ offset2;

phi_mod  = phi_mod + offset3;
phi_mod2 = phi_mod2+ offset4;

phi_mod_gt3 = phi_mod_gt3+ offset5;
phi_mod_gt4 = phi_mod_gt4+ offset6;

diff = phi_mod_gt - phi_prev_mod_gt_load;

%-------------------
%diff = phi_mod - phi_prev_mod;

delta_phi = diff;

% +++++++++++++++++++++++++++++++++++++++++++++ Noise of Phase +++++++++++++++++++++++++++++++++++++++++++++++

% new_sigma = sigma/phi_mod_gt + sigma/phi_mod_gt;
% phi_noise = 0.04*new_sigma.*rand;
% 
% new_sigma2 = sigma/phi_mod_gt2 + sigma/phi_mod_gt2;
% phi_noise2 = 0.04*new_sigma2.*rand;
% 
% new_sigma3 = sigma/phi_mod_gt3 + sigma/phi_mod_gt3;
% phi_noise3 = 0.04*new_sigma3.*rand;
% 
% new_sigma4 = sigma/phi_mod_gt4 + sigma/phi_mod_gt4;
% phi_noise4 = 0.04*new_sigma4.*rand;
% 
% phi_mod_gt  = exp(phi_noise)*phi_mod_gt;
% phi_mod_gt2 = exp(phi_noise2)*phi_mod_gt2;
% % 
% phi_mod_gt3  = exp(phi_noise3)*phi_mod_gt3;
% phi_mod_gt4  = exp(phi_noise4)*phi_mod_gt4;

if isnan(phi_prev_mod_load)
    phi_prev_mod_load = phi_mod2;
end

if isnan(phi_prev_mod_gt_load)
    phi_prev_mod_gt_load = phi_mod_gt;
end

rphase_prev = 3*10^8/(4*pi)*(phi_prev_mod_gt - phi_prev_mod2_gt)/Delta_f;
rphase  = 3*10^8/(2*pi)*1/2*(1*(phi_mod_gt - phi_mod_gt2)+1*(phi_mod_gt4 - phi_mod_gt3))/Delta_f;

rdot   = -lambda/(2*pi)*(phi_mod_gt - phi_prev_mod_gt_load)*1/T;
rdotgt = -lambda/(2*pi)*diff*1/T;

if abs(rdot) > 1
    rdot = 0;
end
if abs(rdotgt) > 1
    rdotgt = 0;
end
if abs(rphase) > 3
    rphase = rphase_prev;
end

end