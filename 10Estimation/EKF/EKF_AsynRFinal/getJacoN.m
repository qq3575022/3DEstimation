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

    
    
%if index == 5
for j = 1 : length(x)
    e = I(:,j);
    E = getEF(y,x-dx*e, factor);
    Edx = getEF(y,x+dx*e, factor);
    J(:,j) = (E-Edx)/(2*dx);
end

%end
    
    

end