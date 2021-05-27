function m2 = getm2(D)

m2 = NaN(1, length(D) - 10000);
%m4 = NaN(1, length(D) - 10000);


for i = 1:1:length(D) - 10000
   % moment1 = mean(D(i:i+10000));
    moment2 = mean(D(i:i+10000).^2);
    moment4 = mean(D(i:i+10000).^4);
   % moment8 = mean(D(i:i+10000).^8);

    m2(i) = (-2*moment2^2 + moment4 - moment2*sqrt(2*moment2^2 - moment4))/(moment2^2 - moment4);
    %(moment2 - moment1.^2)/4 - 1;%abs(moment1^2/(moment2 - moment1^2));
    %m4(i) = abs((-2*moment4^2 + moment8 - moment4*sqrt(abs(2*moment4^2 - moment8)))/(moment4^2 - moment8));
end

end