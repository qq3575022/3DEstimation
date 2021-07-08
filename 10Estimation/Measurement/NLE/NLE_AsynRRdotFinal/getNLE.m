function E = getNLE(y, x, N, index, len, time, m, factor)

HN = NaN(length(y),1);

%factor = 1000;

% length(y)
% len
% index
% m

for i = 1:1:N

     Gain = getF(time, m, m+i-1);

    if index(i) == 1

        HN(len(i)+1:len(i+1)) = getHxkF(Gain*x, factor);

    end

end
EN = NaN(length(y),1);

for j = 1:1:length(y)
    EN(j) = y(j) - HN(j);
end

E = EN.^2;
%EN
end