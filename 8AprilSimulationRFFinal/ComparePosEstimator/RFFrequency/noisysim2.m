function [v, v2, r, r2, r3, rdot, rdot2, phi_mod_gt, phi_mod_gt2, phi_mod, phi_mod2, offset1, offset2, offset3, offset4] = noisysim2(x,f,Gt,M,X,PT,GT,GR,R,sigma,k,z,z_prev,phi_prev_mod_gt_load, phi_prev_mod_load, T, Hmeas, offset1, offset2, offset3, offset4)
% ============================================= Magnitude =================================================

lambda = 3*10^8/f; lambda2 = 3*10^8/(f+0.1*10^6);
xcoord = z(1,k)-x(1);  ycoord = z(2,k)-x(2);   zcoord = z(3,k)-x(3);

% ============================================= Magnitude =================================================
direct        = sqrt((xcoord)^2+(ycoord)^2+(zcoord)^2);
mul           = 1/direct*exp(-i*2*pi*direct/lambda);
mul2          = 1/direct*exp(-i*2*pi*direct/lambda2);

direct_prev   = sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2+(z_prev(3,k)-x(3))^2);
mul_prev      = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda);
mul_prev2     = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda2);

drefl1        = sqrt((x(1) + 3)^2   + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2)      + sqrt((z(1,k) + 3)^2      + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef1       = 1/drefl1*exp(-i*2*pi*drefl1/lambda);
mulRef12      = 1/drefl1*exp(-i*2*pi*drefl1/lambda2);

drefl1_prev   = sqrt((x(1) + 3)^2   + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2) + sqrt((z_prev(1,k) + 3)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef1_prev  = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda);
mulRef12_prev = 1/drefl1_prev*exp(-i*2*pi*drefl1_prev/lambda2);

drefl2        = sqrt((4.32- x(1))^2   + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2)      + sqrt((4.32 - z(1,k))^2 + ((x(2) - z(2,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef2       = 1/drefl2*exp(-i*2*pi*drefl2/lambda);
mulRef22      = 1/drefl2*exp(-i*2*pi*drefl2/lambda2);

drefl2_prev   = sqrt((4.32- x(1))^2   + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((4.32 - z_prev(1,k))^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef2_prev  = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda);
mulRef22_prev = 1/drefl2*exp(-i*2*pi*drefl2_prev/lambda2);

drefl3        = sqrt((x(2) + 3.36)^2   + ((x(1) - z(1,k))/2)^2     + ((x(3) - z(3,k))/2)^2)      + sqrt((z(2,k) + 3.36)^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef3       = 1/drefl3*exp(-i*2*pi*drefl3/lambda);
mulRef32      = 1/drefl3*exp(-i*2*pi*drefl3/lambda2);

drefl3_prev   = sqrt((x(2) + 3.36)^2   + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((z_prev(2,k) + 3.36)^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef3_prev  = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda);
mulRef32_prev = 1/drefl3_prev*exp(-i*2*pi*drefl3_prev/lambda2);

drefl4        = sqrt((3.44 - x(2))^2   + ((x(1) - z(1,k))/2)^2     + ((x(3) - z(3,k))/2)^2)      + sqrt((3.44 - z(2,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
mulRef4       = 1/drefl4*exp(-i*2*pi*drefl4/lambda);
mulRef42      = 1/drefl4*exp(-i*2*pi*drefl4/lambda2);

drefl4_prev   = sqrt((3.44 - x(2))^2   + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2)+ sqrt((3.44 - z_prev(2,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;
mulRef4_prev  = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda);
mulRef42_prev = 1/drefl4_prev*exp(-i*2*pi*drefl4_prev/lambda2);

drefl5        = sqrt((x(3) + 1)^2   + ((x(2) - z(2,k))/2)^2     + ((x(1) - z(1,k))/2)^2)      + sqrt((z(3,k) + 1)^2 + ((x(2) - z(2,k))/2)^2+ ((x(1) - z(1,k))/2)^2);%direct + lambda*50;
mulRef5       = 1/drefl5*exp(-i*2*pi*drefl5/lambda);
mulRef52      = 1/drefl5*exp(-i*2*pi*drefl5/lambda2);

drefl5_prev   = sqrt((x(3) + 1)^2   + ((x(2) - z_prev(2,k))/2)^2+ ((x(1) - z_prev(1,k))/2)^2) + sqrt((z_prev(3,k) + 1)^2 + ((x(2) - z_prev(2,k))/2)^2+ ((x(1) - z_prev(1,k))/2)^2);%direct + lambda*50;
mulRef5_prev  = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda);
mulRef52_prev = 1/drefl5_prev*exp(-i*2*pi*drefl5_prev/lambda2);

drefl6        = sqrt((2 - x(3))^2  + ((x(1) - z(1,k))/2)^2     + ((x(2) - z(2,k))/2)^2)       + sqrt((2 - z(3,k))^2 + ((x(1) - z(1,k))/2)^2+ ((x(2) - z(2,k))/2)^2);%direct + lambda*50;
mulRef6       = 1/drefl6*exp(-i*2*pi*drefl6/lambda);
mulRef62      = 1/drefl6*exp(-i*2*pi*drefl6/lambda2);

drefl6_prev   = sqrt((2 - x(3))^2  + ((x(1) - z_prev(1,k))/2)^2 + ((x(2) - z_prev(2,k))/2)^2) + sqrt((2 - z_prev(3,k))^2 + ((x(1) - z_prev(1,k))/2)^2+ ((x(2) - z_prev(2,k))/2)^2);%direct + lambda*50;
mulRef6_prev  = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda);
mulRef62_prev = 1/drefl6_prev*exp(-i*2*pi*drefl6_prev/lambda2);

H  = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul+ 0.2*mulRef1 + 0.2*mulRef2 + 0.2*mulRef3 + 0.2*mulRef4 + 0.2*mulRef5 +0.2*mulRef6))^4));%  + 0.6*mulSum3 %mul+ 0.25* mulSum + 0.25*mulSum2
H2 = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul))^4));%mul+ 0.25* mulSum + 0.25*mulSum2

% ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++


H = H^2/abs(H);
pd = makedist('Rician','s',sqrt(H),'sigma',0.2*sqrt(H));
H = random(pd).^2;
v = H;%random(pd).^2;
v2 = H2;

unirand = rand;

r  = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*H^2))^(-1/4);
r2 = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*H2^2))^(-1/4);
r3 = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*Hmeas^2))^(-1/4);
%%

% % ============================================== Phase ====================================================
%------------------------------

phi_prev_mod_gt  = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod_gt       = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
%phi_prev_mod_gt = 4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2+(z_prev(3,k)-x(3))^2)/lambda;
%phi_mod_gt      = 4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

phi_prev_mod2_gt = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul_prev2)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod_gt2      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul2)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
%phi_mod_gt2     = 4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

phi_prev_mod = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul_prev + 0.2*mulRef1_prev+ 0.2*mulRef2_prev+ 0.2*mulRef3_prev+ 0.2*mulRef4_prev+ 0.2*mulRef5_prev+ 0.2*mulRef6_prev)^4));
phi_mod      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul      + 0.2*mulRef1     + 0.2*mulRef2     + 0.2*mulRef3     + 0.2*mulRef4     + 0.2*mulRef5     + 0.2*mulRef6)^4));
%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);

phi_prev_mod2= angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul_prev2+ 0.2*mulRef12_prev+ 0.2*mulRef22_prev+ 0.2*mulRef32_prev+ 0.2*mulRef42_prev+ 0.2*mulRef52_prev+ 0.2*mulRef62_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod2     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul2     + 0.2*mulRef12     + 0.2*mulRef22     + 0.2*mulRef32     + 0.2*mulRef42     + 0.2*mulRef52     + 0.2*mulRef62)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

if abs(phi_prev_mod_gt - phi_mod_gt) > 1
    offset1 = offset1 + pi*sign(phi_prev_mod_gt - phi_mod_gt);
end

if abs(phi_prev_mod2_gt - phi_mod_gt2) > 1
    offset2 = offset2 + pi*sign(phi_prev_mod2_gt - phi_mod_gt2);
end

if abs(phi_prev_mod - phi_mod) > 1
    offset3 = offset3 + pi*sign(phi_prev_mod - phi_mod);
end

if abs(phi_prev_mod2 - phi_mod2) > 1
    offset4 = offset4 + pi*sign(phi_prev_mod2 - phi_mod2);
end

phi_mod_gt  = phi_mod_gt + offset1;
phi_mod_gt2 = phi_mod_gt2+ offset2;

phi_mod  = phi_mod + offset3;
phi_mod2 = phi_mod2+ offset4;

diff = phi_mod - phi_prev_mod_load;

%-------------------
%diff = phi_mod - phi_prev_mod;

delta_phi = diff;

% +++++++++++++++++++++++++++++++++++++++++++++ Noise of Phase +++++++++++++++++++++++++++++++++++++++++++++++

new_sigma = sigma/phi_mod + sigma/phi_prev_mod;
phi_noise = 8000*new_sigma.*rand;

phi_mod_gt = exp(phi_noise)*phi_mod_gt;
phi_mod_gt2 = exp(phi_noise)*phi_mod_gt2;
% phi_mod = phi_mod
% phi_prev_mod_load =  phi_prev_mod_load

% rdot2 = lambda/(4*pi)*delta_phi2*1/T;
rdot  = 3*10^8/(4*pi)*(phi_mod -phi_mod2)/(0.1*10^6);%lambda/(4*pi)*phi_mod*1/T;

%rdot  = 2*3*10^8/(4*pi)*(phi_mod -phi_mod2)/(0.1);%lambda/(4*pi)*phi_mod*1/T;
rdot2 = lambda/(4*pi)*diff*1/T;

if abs(rdot2) > 2
    rdot2 = 0;
end
if abs(rdot) > 3
    rdot = 0;
end
end
%---------------------------------
% phi_prev_mod = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul_prev))^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
% 
% phi_mod      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul + 0.2*mulRef1))^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% 
% %phi_mod2     = angle(sqrt(2*R*PT*GT*GR*Gt^2*lambda2^2*X^2*M)/(4*pi)*(mul2 + 0.2*mulRef12));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% phi_mod2      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul2 + 0.2*mulRef12))^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% 
% %phi_mod_gt      = angle(sqrt(2*R*PT*GT*GR*Gt^2*lambda^2*X^2*M)/(4*pi)*(mul))%  + 0.2*mulRef1));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% phi_mod2      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul))^4));
% 
% %phi_mod2_gt     = angle(sqrt(2*R*PT*GT*GR*Gt^2*lambda2^2*X^2*M)/(4*pi)*(mul2))%+ 0.2*mulRef12));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
% phi_mod2_gt      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/((4*pi)*(mul2))^4));

%----------------------------------

% global l; global l1; global l2;global l3;
% 
% if index == 0
%     if (phi_prev_mod - phi_mod ) > 4.5
%         l = l + 1;
%         phi_prev_mu = phi_prev_mod + (l-1)*2*pi;
%         phi_mu = phi_mod + l*2*pi;
%     elseif phi_mod > 5 && phi_prev_mod < 1 
%         l = l - 1;
%         phi_prev_mu = phi_prev_mod + l*2*pi;
%         phi_mu = phi_mod + (l-1)*2*pi;
%     else
%         phi_prev_mu = phi_prev_mod + l*2*pi;
%         phi_mu = phi_mod + l*2*pi;
%     end
%     
% elseif index == 1
%     if (phi_prev_mod - phi_mod)> 1.8 %5 && phi_mod < 1 
%         l1 = l1 + 1;
%         phi_prev_mu = phi_prev_mod + (l1-1)*2*pi;
%         phi_mu = phi_mod + l1*2*pi;
% 
%     elseif (phi_mod - phi_prev_mod)> 2.8 %phi_prev_mod < 1 && phi_mod > 5
%         l1 = l1 - 1;
%         phi_prev_mu = phi_prev_mod + l1*2*pi;
%         phi_mu = phi_mod + (l1-1)*2*pi;
% 
%     else
%         phi_prev_mu = phi_prev_mod + l1*2*pi;
%         phi_mu = phi_mod + l1*2*pi;
%     end
%     
% elseif index == 2
%     if k < 0.47*length(z)
%         if (phi_prev_mod - phi_mod)> 2.9 %5 && phi_mod < 1 
%             l2 = l2 + 1;
%             phi_prev_mu = phi_prev_mod + (l2-1)*2*pi;
%             phi_mu = phi_mod + l2*2*pi;
% 
%         else
%             phi_prev_mu = phi_prev_mod + l2*2*pi;
%             phi_mu = phi_mod + l2*2*pi;
%         end
%     else
%         
%         if (phi_mod - phi_prev_mod)> 3.0 %phi_prev_mod < 1 && phi_mod > 5
%             l2 = l2 - 1;
%             phi_prev_mu = phi_prev_mod + l2*2*pi;
%             phi_mu = phi_mod + (l2-1)*2*pi;
% 
%         else
%             phi_prev_mu = phi_prev_mod + l2*2*pi;
%             phi_mu = phi_mod + l2*2*pi;
%         end
%     end
%     
% else 
%     if k < 0.5*length(z)
%         if (phi_mod - phi_prev_mod)> 2.3 %phi_prev_mod < 1 && phi_mod > 5
%             l3 = l3 - 1;
%             phi_prev_mu = phi_prev_mod + l3*2*pi;
%             phi_mu = phi_mod + (l3-1)*2*pi;
%             
%         elseif (phi_prev_mod - phi_mod)> 5.0 %5 && phi_mod < 1 
%             l3 = l3 + 1;
%             phi_prev_mu = phi_prev_mod + (l3-1)*2*pi;
%             phi_mu = phi_mod + l3*2*pi;
%         else
%             phi_prev_mu = phi_prev_mod + l3*2*pi;
%             phi_mu = phi_mod + l3*2*pi;
%         end
%         
%     else
%         if (phi_prev_mod - phi_mod)> 3.0%5 && phi_mod < 1 
%             l3 = l3 + 1;
%             phi_prev_mu = phi_prev_mod + (l3-1)*2*pi;
%             phi_mu = phi_mod + l3*2*pi;
%    
%         else
%             phi_prev_mu = phi_prev_mod + l3*2*pi;
%             phi_mu = phi_mod + l3*2*pi;
%         end
%                
%     end
% end

%phi_prev_conc = phi_prev_mu;  %icdf('Normal',uniphase,phi_prev_mu,sigma);
%phi_conc      = phi_mu;       %icdf('Normal',uniphase2,phi_mu,sigma);