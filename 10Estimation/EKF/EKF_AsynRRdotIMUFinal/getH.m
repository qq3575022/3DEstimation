function H = getH(x, index, factor)

%H = getHxkA(x);
if index == 1
    H = getHxkFGA(x, factor);

elseif index == 2
    H = getHxkFG(x, factor);

elseif index == 3
    H = getHxkFA(x, factor); 
    
elseif index == 4
    H = getHxkGA(x, factor);
    
elseif index == 5
    H = getHxkF(x, factor);

elseif index == 6
    H = getHxkA(x, factor); 
    
elseif index == 7    
    H = getHxkG(x, factor); 
end

end