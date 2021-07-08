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
        HN(len(i)+1:len(i+1)) = getHxkFGA(Gain*x, factor);

    elseif index(i) == 2
        HN(len(i)+1:len(i+1)) = getHxkFG(Gain*x, factor);

    elseif index(i) == 3
        HN(len(i)+1:len(i+1)) = getHxkFA(Gain*x, factor);

    elseif index(i) == 4
        HN(len(i)+1:len(i+1)) = getHxkGA(Gain*x, factor);

    elseif index(i) == 5
        HN(len(i)+1:len(i+1)) = getHxkF(Gain*x, factor);

    elseif index(i) == 6
        HN(len(i)+1:len(i+1)) = getHxkA(Gain*x, factor);

    elseif index(i) == 7
        HN(len(i)+1:len(i+1)) = getHxkG(Gain*x, factor);
    end

end
EN = NaN(length(y),1);

for j = 1:1:length(y)
    EN(j) = y(j) - HN(j);
end

E = EN.^2;
%EN
end