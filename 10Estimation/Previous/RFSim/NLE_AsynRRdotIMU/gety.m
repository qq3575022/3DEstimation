function [y, ii, jj, kk, index, len] = gety(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m, N)
    
    y = []; length = 0; len = [0]; index = [];
    
    for l = 1:1:N
        
        [y1, i, j, k, index1, len1] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m+l-1);
        
        if l == 1
            ii = i; jj = j; kk = k;
        end
        y = [y;y1];
        
        %index1
        
        length = length + len1;
        len = [len; length];
        
        index = [index; index1];
%         ii = i; jj = j; kk = k; 
%         [y2, i, j, k, index2, len2] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m+1);
%         [y3, i, j, k, index3, len3] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m+2);
%         [y4, i, j, k, index4, len4] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m+3);
%         [y5, i, j, k, index5, len5] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,mag_data_x, gyro_data_x, mag_data_y, gyro_data_y, mag_data_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, mag_time, time, i, j, k, m+4);
    end
  
    
%     y = [y1;y2;y3;y4;y5];
%     index = [index1, index2, index3, index4, index5];
%     len = [0, len1, len1 + len2, len1 + len2 + len3, len1 + len2 + len3 + len4, len1 + len2 + len3 + len4 + len5];
end