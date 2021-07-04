function [y, l, index1, len1] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, RFtime, l, m)
% m = m
% i = i
% j = j
% l = l
% timeTt = time(m)
% accTt = acc_time(i)
% gyroTt = gyro_time(j)
% RF_time = RFtime(l)

factor = 1;
  
if  l <= length(RFtime) 
  index = 1;
  len = 8;
  y = NaN(8,1);
  y(1:8) = [r_sim1(l), rdot_sim1(l), r_sim2(l), rdot_sim2(l), r_sim3(l), rdot_sim3(l), r_sim4(l), rdot_sim4(l)];
  l = l + 1;
  

else
    y = nan;
    len = -1;
    index = -1;
  
end

len1 = len;
index1 = index;


end