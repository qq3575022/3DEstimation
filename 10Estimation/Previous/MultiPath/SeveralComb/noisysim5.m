function [v, v2, r, r2, rdot, rdot2] = noisysim5(x,f,Gt,M,X,PT,GT,GR,R,sigma,index,k,z,z_prev,T,time1,multi)

lambda1  = 3*10^8/f;

% ============================================= Magnitude =================================================

direct      = sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2);
direct_prev = sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2+(z_prev(3,k)-x(3))^2);

mul1        = 1/direct*exp(-i*2*pi*direct/lambda1);
mul_prev    = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda1);


% multi path
d3          = sqrt((x(1) + 2)^2 + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2)      + sqrt((z(1,k) + 2)^2      + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
d3_prev     = sqrt((x(1) + 2)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((z_prev(1,k) + 2)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;

d4          = sqrt((3- x(1))^2  + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2)      + sqrt((3 - z(1,k))^2      + ((x(2) - z(2,k))/2)^2      + ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
d4_prev     = sqrt((3- x(1))^2  + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((3 - z_prev(1,k))^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;

d5          = sqrt((x(2) + 2)^2 + ((x(1) - z(1,k))/2)^2      + ((x(3) - z(3,k))/2)^2)      + sqrt((z(2,k) + 2)^2      + ((x(1) - z(1,k))/2)^2      + ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
d5_prev     = sqrt((x(2) + 2)^2 + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((z_prev(2,k) + 2)^2 + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;

d6          = sqrt((3 - x(2))^2 + ((x(1) - z(1,k))/2)^2      + ((x(3) - z(3,k))/2)^2)      + sqrt((3 - z(2,k))^2      + ((x(1) - z(1,k))/2)^2      + ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
d6_prev     = sqrt((3 - x(2))^2 + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((3 - z_prev(2,k))^2 + ((x(1) - z_prev(1,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;

d7          = sqrt((x(3) + 1)^2 + ((x(2) - z(2,k))/2)^2      + ((x(1) - z(1,k))/2)^2)      + sqrt((z(3,k) + 1)^2      + ((x(2) - z(2,k))/2)^2      + ((x(1) - z(1,k))/2)^2);%direct + lambda*50;
d7_prev     = sqrt((x(3) + 1)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(1) - z_prev(1,k))/2)^2) + sqrt((z_prev(3,k) + 1)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(1) - z_prev(1,k))/2)^2);%direct + lambda*50;

d8         = sqrt((5 - x(3))^2 + ((x(1) - z(1,k))/2)^2       + ((x(2) - z(2,k))/2)^2)      + sqrt((5 - z(3,k))^2      + ((x(1) - z(1,k))/2)^2      + ((x(2) - z(2,k))/2)^2);%direct + lambda*50;
d8_prev    = sqrt((5 - x(3))^2 + ((x(1) - z_prev(1,k))/2)^2  + ((x(2) - z_prev(2,k))/2)^2) + sqrt((5 - z_prev(3,k))^2 + ((x(1) - z_prev(1,k))/2)^2 + ((x(2) - z_prev(2,k))/2)^2);%direct + lambda*50;



mulSum3      = 1/d3*exp(-i*2*pi*d3/lambda1);  % frequency 1
mulSum3_prev = 1/d3_prev*exp(-i*2*pi*d3_prev/lambda1);

mulSum4      = 1/d4*exp(-i*2*pi*d4/lambda1);
mulSum4_prev = 1/d4_prev*exp(-i*2*pi*d4_prev/lambda1);

mulSum5      = 1/d5*exp(-i*2*pi*d5/lambda1);
mulSum5_prev = 1/d5_prev*exp(-i*2*pi*d5_prev/lambda1);

mulSum6      = 1/d6*exp(-i*2*pi*d6/lambda1);
mulSum6_prev = 1/d6_prev*exp(-i*2*pi*d6_prev/lambda1);

mulSum7      = 1/d7*exp(-i*2*pi*d7/lambda1);
mulSum7_prev = 1/d7_prev*exp(-i*2*pi*d7_prev/lambda1);

mulSum8      = 1/d8*exp(-i*2*pi*d8/lambda1);
mulSum8_prev = 1/d8_prev*exp(-i*2*pi*d8_prev/lambda1);

H      = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda1^2*X^2*M)/(4*pi)^2*(mul1      + 0.2*mulSum3 + 0.2*mulSum4 + 0.2*mulSum5 + 0.1*mulSum6 + 0.1*mulSum7 +0.2*mulSum8)^2));%mul + mulSum  % + 0.25* mulSum3 %mul + mulSum2
H_prev = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda1^2*X^2*M)/(4*pi)^2*(mul_prev  + 0.2*mulSum3_prev + 0.2*mulSum4_prev + 0.2*mulSum5_prev + 0.1*mulSum6_prev + 0.1*mulSum7_prev +0.2*mulSum8_prev)^2));

% ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++
%pd = makedist('Rician','s',sqrt(H),'sigma',0.01*sqrt(H));
%v = random(pd).^2;

v  = H;
v2 = H;

unirand = rand;

r = ((2*PT*GT*GR*Gt^2*lambda1^2*X^2*M*R)/((4*pi)^2*(H^2)))^(1/2);
r2 = r;

% % ============================================== Phase ====================================================

phi_mod1     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda1^2*X^2*M)/(4*pi)^2*(mul1     + 0.2*mulSum3 + 0.2*mulSum4 + 0.2*mulSum5 + 0.1*mulSum6 + 0.1*mulSum7 +0.2*mulSum8)^2));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
phi_prev_mod = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda1^2*X^2*M)/(4*pi)^2*(mul_prev + 0.2*mulSum3_prev + 0.2*mulSum4_prev + 0.2*mulSum5_prev + 0.1*mulSum6_prev + 0.1*mulSum7_prev +0.2*mulSum8_prev)^2));

diff = phi_mod1 - phi_prev_mod;

delta_phi = diff;%phi_conc - phi_prev_conc;

% +++++++++++++++++++++++++++++++++++++++++++++ Noise of Phase +++++++++++++++++++++++++++++++++++++++++++++++

new_sigma = sigma/phi_mod1 + sigma/phi_prev_mod;
phi_noise = 800000*new_sigma.*rand;

delta_phi2 = exp(phi_noise)*delta_phi;


% rdot  = lambda/(4*pi)*diff*1/T;
% rdot2 = lambda/(4*pi)*delta_phi2*1/T;

% rdot  = 3*10^8/(2*pi)*(phi_mod1 -phi_mod2)/(0.1);%lambda/(4*pi)*phi_mod*1/T;

rdot = -lambda1/(4*pi)*diff*1/T;
rdot2 = -lambda1/(4*pi)*diff*1/T;

if abs(rdot) > 1
    rdot = 0;
end

if abs(rdot2) > 1 
    rdot2 = 0;
end

%v2 = icdf('Normal',unirand,v,sigma); %sigma*randn + mu;
%r2 = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*v2^2))^(1/4);

end


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