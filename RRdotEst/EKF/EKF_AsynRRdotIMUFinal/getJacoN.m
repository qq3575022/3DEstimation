function J = getJacoN(y,x,index, factor)

dx = 1e-12;
I = eye(length(x));
J = [ ];
%index
% for j = 1 : length(x)
%         e = I(:,j);
%         E = getEA(y,x-dx*e);
%         Edx = getEA(y,x+dx*e);
%         J(:,j) = (E-Edx)/(2*dx);
% end

    
if index == 1
    for j = 1 : length(x)
        e = I(:,j);
        E = getEFGA(y,x-dx*e, factor);
        Edx = getEFGA(y,x+dx*e, factor);
        J(:,j) = (E-Edx)/(2*dx);
    end

elseif index == 2
    for j = 1 : length(x)
        e = I(:,j);
        E = getEFG(y,x-dx*e, factor);
        Edx = getEFG(y,x+dx*e, factor);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
elseif index == 3
    for j = 1 : length(x)
        e = I(:,j);
        E = getEFA(y,x-dx*e, factor);
        Edx = getEFA(y,x+dx*e, factor);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
elseif index == 4
    for j = 1 : length(x)
        e = I(:,j);
        E = getEGA(y,x-dx*e, factor);
        Edx = getEGA(y,x+dx*e, factor);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
elseif index == 5
    for j = 1 : length(x)
        e = I(:,j);
        E = getEF(y,x-dx*e, factor);
        Edx = getEF(y,x+dx*e, factor);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
elseif index == 6
    for j = 1 : length(x)
        e = I(:,j);
        E = getEA(y,x-dx*e, factor);
        Edx = getEA(y,x+dx*e, factor);
        J(:,j) = (E-Edx)/(2*dx);
    end

elseif index == 7
    for j = 1 : length(x)
        e = I(:,j);
        E = getEG(y,x-dx*e, factor);
        Edx = getEG(y,x+dx*e, factor);
        J(:,j) = (E-Edx)/(2*dx);
    end
    
    

end