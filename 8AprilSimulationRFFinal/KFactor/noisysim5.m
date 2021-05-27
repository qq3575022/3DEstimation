function [v, v2, r, r2, rdot, rdot2] = noisysim5(x,f,Gt,M,X,PT,GT,GR,R,sigma,index,k,z,z_prev,T,time1,multi)

lambda1  = 3*10^8/f;
lambda2  = 3*10^8/(f + 0.1);

% ============================================= Magnitude =================================================

direct      = sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2);
mul1         = 1/direct*exp(-i*2*pi*direct/lambda1);
mul2         = 1/direct*exp(-i*2*pi*direct/lambda2);

direct_prev = sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2+(z_prev(3,k)-x(3))^2);
mul_prev    = 1/direct_prev*exp(-i*2*pi*direct_prev/lambda1);



H      = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda1^2*X^2*M)/(4*pi)^2*(mul1)^2));%mul + mulSum  % + 0.25* mulSum3 %mul + mulSum2
H_prev = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda1^2*X^2*M)/(4*pi)^2*(mul_prev)^2));

% ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++
%pd = makedist('Rician','s',sqrt(H),'sigma',0.01*sqrt(H));
%v = random(pd).^2;

v = H;
v2 = H;

unirand = rand;

r = ((2*PT*GT*GR*Gt^2*lambda1^2*X^2*M*R)/((4*pi)^2*(H^2)))^(1/2);
r2 = r;

% % ============================================== Phase ====================================================

phi_prev_mod = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda1^2*X^2*M)/(4*pi)^2*(mul_prev)^2));
phi_mod1     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda1^2*X^2*M)/(4*pi)^2*(mul1)^2));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
phi_mod2     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^2*X^2*M)/(4*pi)^2*(mul2)^2));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);

diff = phi_mod1 - phi_prev_mod;

delta_phi = diff;%phi_conc - phi_prev_conc;

% +++++++++++++++++++++++++++++++++++++++++++++ Noise of Phase +++++++++++++++++++++++++++++++++++++++++++++++

new_sigma = sigma/phi_mod1 + sigma/phi_prev_mod;
phi_noise = 800000*new_sigma.*rand;

delta_phi2 = exp(phi_noise)*delta_phi;


% rdot  = lambda/(4*pi)*diff*1/T;
% rdot2 = lambda/(4*pi)*delta_phi2*1/T;

rdot  = 3*10^8/(2*pi)*(phi_mod1 -phi_mod2)/(0.1);%lambda/(4*pi)*phi_mod*1/T;


rdot2 = lambda1/(4*pi)*diff*1/T;

if abs(rdot2) > 3
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