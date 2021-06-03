function E = getMulPath(z1, xx, gamma, f, phi)

%phi1 = mod(phi, 2*pi);

y = getPhase(z1, xx, gamma, f);

%y1= mod(y, 2*pi);

E = (y - phi).^2;

end

    
%     if index(i) == 1
%         HN(len(i)+1:len(i+1)) = getHxkMGA(Gain*x);
% 
%     elseif index(i) == 2
%         HN(len(i)+1:len(i+1)) = getHxkGA(Gain*x);
% 
%     elseif index(i) == 3
%         HN(len(i)+1:len(i+1)) = getHxkMG(Gain*x);
% 
%     elseif index(i) == 4
%         HN(len(i)+1:len(i+1)) = getHxkMA(Gain*x);
% 
%     elseif index(i) == 5
%         HN(len(i)+1:len(i+1)) = getHxkA(Gain*x);
% 
%     elseif index(i) == 6
%         HN(len(i)+1:len(i+1)) = getHxkG(Gain*x);
% 
%     elseif index(i) == 7
%         HN(len(i)+1:len(i+1)) = getHxkM(Gain*x);
%     end