function [v, v2, r, r2, rdot, rdot2, diff] = noisysim3(x,f,Gt,M,X,PT,GT,GR,R,sigma,index,k,z,z_prev,T,time1,multi)

lambda  = 3*10^8/f;
lambda2 = 3*10^8/(f+0.1);
% ============================================= Magnitude =================================================
direct       = 1/3*(sqrt((z(1,k)-x(1))^2 + (z(2,k)-x(2))^2 + (z(3,k)-x(3))^2) + sqrt((z(1,k-1)-x(1))^2 + (z(2,k-1)-x(2))^2 + (z(3,k-1)-x(3))^2) +sqrt((z(1,k-2)-x(1))^2 + (z(2,k-2)-x(2))^2 + (z(3,k-2)-x(3))^2));
direct_prev  = 1/3*(sqrt((z_prev(1,k)-x(1))^2 + (z_prev(2,k)-x(2))^2 + (z_prev(3,k)-x(3))^2) + sqrt((z_prev(1,k-1)-x(1))^2 + (z_prev(2,k-1)-x(2))^2 + (z_prev(3,k-1)-x(3))^2) + sqrt((z_prev(1,k-2)-x(1))^2 + (z_prev(2,k-2)-x(2))^2 + (z_prev(3,k-2)-x(3))^2));

mul          =  1/direct*exp(-i*2*pi*direct/lambda);  % frequency 1
mul2         =  1/direct*exp(-i*2*pi*direct/lambda2); % frequency 2
mul_prev     =  1/direct_prev*exp(-i*2*pi*direct_prev/lambda);

% multipath 

d3           = sqrt((x(1) + 2)^2 + ((x(2) - z(2,k))/2)^2 + ((x(3) - z(3,k))/2)^2)+ sqrt((z(1,k) + 2)^2 + ((x(2) - z(2,k))/2)^2  + ((x(3) - z(3,k))/2)^2);%direct + lambda*50;
d3_prev      = sqrt((x(1) + 2)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2) + sqrt((z_prev(1,k) + 2)^2 + ((x(2) - z_prev(2,k))/2)^2 + ((x(3) - z_prev(3,k))/2)^2);%direct + lambda*50;

mulSum3      = 1/d3*exp(-i*2*pi*d3/lambda);  % frequency 1
mulSum4      = 1/d3*exp(-i*2*pi*d3/lambda2); % frequency 2
mulSum3_prev = 1/d3_prev*exp(-i*2*pi*d3_prev/lambda);

H = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul)^4));%mul + mulSum  % + 0.25* mulSum3 %mul + mulSum2
H = H./sqrt(H);

% ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++
H2 = abs(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul)^4));%mul+ 0.25* mulSum + 0.25*mulSum2
H2 = H2./sqrt(H2);
% ++++++++++++++++++++++++++++++++++++++++++ Noise of Magnitude +++++++++++++++++++++++++++++++++++++++++++
pd = makedist('Rician','s',sqrt(H),'sigma',0.01*sqrt(H));
 
%v = random(pd).^2;
v = random(pd).^2;
v2 = H2;
unirand = rand;


%r = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*(H^2)))^(1/4);
r2 = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*(H^2)))^(1/4);

% % ============================================== Phase ====================================================

phi_prev_mod = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul_prev + 0.2*mulSum3_prev)^4));%mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod      = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi)^4*(mul      + 0.2*mulSum3)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda; %4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
phi_mod2     = angle(sqrt((2*R*PT*GT*GR*Gt^2*lambda2^4*X^2*M)/(4*pi)^4*(mul2    + 0.2*mulSum4)^4));%4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2+(z(3,k)-x(3))^2)/lambda2;%mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);


diff = phi_mod - phi_prev_mod;

delta_phi = diff;%phi_conc - phi_prev_conc;

% +++++++++++++++++++++++++++++++++++++++++++++ Noise of Phase +++++++++++++++++++++++++++++++++++++++++++++++

new_sigma = sigma/phi_mod + sigma/phi_prev_mod;
phi_noise = 800000*new_sigma.*rand;

delta_phi2 = exp(phi_noise)*delta_phi;


% rdot  = lambda/(4*pi)*diff*1/T;
% rdot2 = lambda/(4*pi)*delta_phi2*1/T;

r  = 3*10^8/(4*pi)*(phi_mod -phi_mod2)/(0.1);%lambda/(4*pi)*phi_mod*1/T;
rdot = -lambda/(4*pi)*diff*1/T;

rdot2 = rdot;

if abs(rdot) > 3
    rdot = 0;
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