function [y, len] = gety(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, m)
    
    [y1, len1] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,m);
    %ii = i; jj = j; kk = k; 
    [y2, len2] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, m+1);
    [y3, len3] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,m+2);
    [y4, len4] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,m+3);
    [y5, len5] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, m+4);

    y = [y1;y2;y3;y4;y5];
    %index = [index1, index2, index3, index4, index5];
    len = [0, len1, len1 + len2, len1 + len2 + len3, len1 + len2 + len3 + len4, len1 + len2 + len3 + len4 + len5];
end