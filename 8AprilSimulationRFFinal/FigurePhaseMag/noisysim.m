function [v, v2, r, r2, rdot, rdot2, phase, phi_mod_gt, phi_mod2_gt, phi_mod, phi_mod2, offset1, offset2, offset3, offset4] = noisysim(x,f,Gt,X,PT,GT,R,sigma,rice,k,z,z_prev,phi_prev_mod_gt_load, phi_prev_mod_load,T,offset1, offset2, offset3, offset4)
% ============================================= Magnitude =================================================

Deltaf = 10^6; lambda = 3*10^8/f; lambda2 = 3*10^8/(f+Deltaf); lambda3 = 3*10^8/(f-Deltaf);
xcoord = z(1,k)-x(1);  ycoord = z(2,k)-x(2);   zcoord = z(3,k)-x(3);

% ============================================= Magnitude =================================================
direct = sqrt((xcoord)^2+(ycoord)^2+(zcoord)^2);
mul    = 1/direct*exp(-i*2*pi*direct/lambda);
mul2   = 1/direct*exp(-i*2*pi*direct/lambda2);

direct_prev =  sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2+(z_prev(3,k)-x(3))^2);
mul_prev    = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda);
mul_prev2   = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda2);
mul_prev3   = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda3);


drefl1 = sqrt((x(1) + 2)^2 + ((x(2) - z(2,k))/2)^2 + ((x(3) - z(3,k))/2)^2) + sqrt((z(1,k) + 2)^2 + ((x(2) - z(2,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef1  = 1/drefl1*exp(-i*2*pi*drefl1/lambda);
mulRef12 = 1/drefl1*exp(-i*2*pi*drefl1/lambda2);

drefl1_prev = sqrt((x(1) + 2)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((z_prev(1,k) + 2)^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef1_prev  = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda);
mulRef12_prev = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda2);

drefl2 = sqrt((3- x(1))^2 + ((x(2) - z(2,k))/2)^2 + ((x(3) - z(3,k))/2)^2) + sqrt((3 - z(1,k))^2 + ((x(2) - z(2,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef2 = 1/drefl2*exp(-i*2*pi*drefl2/lambda);
mulRef22 = 1/drefl2*exp(-i*2*pi*drefl2/lambda2);

drefl2_prev = sqrt((3- x(1))^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((3 - z_prev(1,k))^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef2_prev = 1/drefl2_prev*exp(-i*2*pi*drefl2_prev/lambda);
mulRef22_prev = 1/drefl2_prev*exp(-i*2*pi*drefl2_prev/lambda2);

drefl3 = sqrt((x(2) + 2)^2 + ((x(1) - z(1,k))/2)^2 + ((x(3) - z(3,k))/2)^2) + sqrt((z(2,k) + 2)^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef3 = 1/drefl3*exp(-i*2*pi*drefl3/lambda);
mulRef32 = 1/drefl3*exp(-i*2*pi*drefl3/lambda2);

drefl3_prev = sqrt((x(2) + 2)^2 + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((z_prev(2,k) + 2)^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef3_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda);
mulRef32_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda2);

drefl4 = sqrt((3 - x(2))^2 + ((x(1) - z(1,k))/2)^2 + ((x(3) - z(3,k))/2)^2) + sqrt((3 - z(2,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef4 = 1/drefl4*exp(-i*2*pi*drefl4/lambda);
mulRef42 = 1/drefl4*exp(-i*2*pi*drefl4/lambda2);

drefl4_prev = sqrt((3 - x(2))^2 + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((3 - z_prev(2,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef4_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda);
mulRef42_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda2);

drefl5 = sqrt((x(3) + 1)^2 + ((x(2) - z(2,k))/2)^2 + ((x(1) - z(1,k))/2)^2) + sqrt((z(3,k) + 1)^2 + ((x(2) - z(2,k))/2)^2+ ((x(1) - z(1,k))/2)^2);%direct + lambda*50;
mulRef5 = 1/drefl5*exp(-i*2*pi*drefl5/lambda);
mulRef52 = 1/drefl5*exp(-i*2*pi*drefl5/lambda2);

drefl5_prev = sqrt((x(3) + 1)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(1) - z_prev(1,k))/2)^2) + sqrt((z_prev(3,k) + 1)^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(1) - z_prev(1,k))/2)^2);%direct + lambda*50;
mulRef5_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda);
mulRef52_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda2);

drefl6 = sqrt((5 - x(3))^2 + ((x(1) - z(1,k))/2)^2 + ((x(2) - z(2,k))/2)^2) + sqrt((5 - z(3,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(2) - z(2,k))/2)^2);%direct + lambda*50;
mulRef6 = 1/drefl6*exp(-i*2*pi*drefl6/lambda);
mulRef62 = 1/drefl6*exp(-i*2*pi*drefl6/lambda2);

drefl6_prev = sqrt((5 - x(3))^2 + ((x(1) - z_prev(1,k))/2)^2 + ((x(2) - z_prev(2,k))/2)^2) + sqrt((5 - z_prev(3,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(2) - z_prev(2,k))/2)^2);%direct + lambda*50;
mulRef6_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda);
mulRef62_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda2);

H  = abs(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul+ 0.1*mulRef1 + 0.1*mulRef2 + 0.1*mulRef3 + 0.1*mulRef4 + 0.1*mulRef5 +0.1*mulRef6));%  + 0.6*mulSum3 %mul+ 0.25* mulSum + 0.25*mulSum2
H2 = abs(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul));%mul+ 0.25* mulSum + 0.25*mulSum2

% ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++

pd = makedist('Rician','s',sqrt(H),'sigma',rice*sqrt(H));
 
v = random(pd).^2;
v2 = H2;

unirand = rand;

r  = sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi*v);
r2 = sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi*H2);
%%

% % ============================================== Phase ====================================================

phi_prev_mod  = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul_prev   + 0.1*mulRef1_prev + 0.1*mulRef2_prev + 0.1*mulRef3_prev + 0.1*mulRef4_prev + 0.1*mulRef5_prev + 0.1*mulRef6_prev));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_prev_mod2 = angle(sqrt(2*R*PT*GT*Gt*lambda2^2*X)/(4*pi)*(mul_prev2 + 0.1*mulRef12_prev+0.1*mulRef22_prev +0.1*mulRef32_prev+0.1*mulRef42_prev+0.1*mulRef52_prev+0.1*mulRef62_prev));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);

phi_prev_mod_gt  = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul_prev));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_prev_mod2_gt = angle(sqrt(2*R*PT*GT*Gt*lambda2^2*X)/(4*pi)*(mul_prev2));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);


phi_prev_mod3 = angle(sqrt(2*R*PT*GT*Gt*lambda3^2*X)/(4*pi)*(mul_prev3));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);

phi_mod      = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul + 0.1*mulRef1 + 0.1*mulRef2 + 0.1*mulRef3+0.1*mulRef4+0.1*mulRef5+0.1*mulRef6));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
phi_mod2     = angle(sqrt(2*R*PT*GT*Gt*lambda2^2*X)/(4*pi)*(mul2 + 0.1*mulRef12+0.1*mulRef22+0.1*mulRef32+0.1*mulRef42+0.1*mulRef52+0.1*mulRef62));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
phase = phi_mod;
% Angle
phi_mod_gt      = angle(sqrt(2*R*PT*GT*Gt*lambda^2*X)/(4*pi)*(mul));%  + 0.2*mulRef1));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
phi_mod2_gt     = angle(sqrt(2*R*PT*GT*Gt*lambda2^2*X)/(4*pi)*(mul2));%+ 0.2*mulRef12));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

%phi_prev_mod    = 4*pi*direct_prev/lambda;
%phi_mod_gt      = 4*pi*direct/lambda;%  + 0.2*mulRef1));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
%phi_mod2_gt     = 4*pi*direct/lambda2;%+ 0.2*mulRef12));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

if abs(phi_prev_mod_gt - phi_mod_gt) > 2
    offset1 = offset1 + phi_prev_mod_gt - phi_mod_gt;
end

if abs(phi_prev_mod2_gt - phi_mod2_gt) > 2
    offset2 = offset2 + phi_prev_mod2_gt - phi_mod2_gt;
end

if abs(phi_prev_mod - phi_mod) > 2
    offset3 = offset3 + phi_prev_mod - phi_mod;
end

if abs(phi_prev_mod2 - phi_mod2) > 2
    offset4 = offset4 + phi_prev_mod2 - phi_mod2;
end

phi_mod_gt  = phi_mod_gt + offset1;
phi_mod2_gt = phi_mod2_gt+ offset2;

phi_mod  = phi_mod + offset3;
phi_mod2 = phi_mod2+ offset4;

diff = phi_mod - phi_prev_mod_load;

delta_phi = diff;

% +++++++++++++++++++++++++++++++++++++++++++++ Noise of Phase +++++++++++++++++++++++++++++++++++++++++++++++

new_sigma = sigma/phi_mod + sigma/phi_prev_mod;
phi_noise = 800000*new_sigma.*rand;

delta_phi2 = exp(phi_noise)*delta_phi;


% rdot2 = lambda/(4*pi)*delta_phi2*1/T;

% lambda
% lambda2
% direct
% phi_mod_gt- phi_mod2_gt

rdot  = 3*10^8/(2*pi)*(phi_prev_mod_load -phi_mod2)/(Deltaf);%lambda/(4*pi)*phi_mod*1/T;
rdot2 = 3*10^8/(2*pi)*(phi_mod -phi_mod2)/(Deltaf);%lambda/(4*pi)*diff*1/T;

if abs(rdot2) > 10
    rdot2 = 0;
end

if abs(rdot) > 10
    rdot = 0;
end

end