clc, clear, close all

data=readtable('2.csv','Delimiter', ',');  g=9.7953;
load('acc_t.mat'); load('std_static.mat');


sensor_acc_check  = ismember(data.Var2,'ACC_UN');   find_acc_data     = find(sensor_acc_check == 1); 

data_x = table2array(data(:,3));  data_y = table2array(data(:,5));  data_z = table2array(data(:,4));


acc_data_x = data_x(find_acc_data);acc_data_y = data_y(find_acc_data);acc_data_z = data_z(find_acc_data);

figure
plot(acc_t, std_static); hold on; plot(acc_t,acc_data_x); hold on; plot(acc_t,acc_data_y);hold on; plot(acc_t,acc_data_z);