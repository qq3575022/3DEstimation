function [phi,phi2, offset1, offset2] = concatenate(phi1, phi21, k, offset1, offset2)

if abs(phi1(k-1) - phi1(k)) > 1
    k
    phi1(k-1)
    phi1(k)
    offset1 = offset1 + phi1(k-1) - phi1(k)
%     if phi1(k-1) > phi1(k)
%         phi = phi1(k) + pi;
%     else
%         phi = phi1(k) - pi;
%     end
end

if abs(phi21(k-1) - phi21(k)) > 1
    k
    phi21(k-1)
    phi21(k)
    offset2 = offset2 + phi21(k-1) - phi21(k)
    
%     if phi21(k-1) > phi21(k)
%         phi2 = phi21(k) + pi;
%     else
%         phi2 = phi21(k) - pi;
%     end
end

phi  = phi1(k) + offset1;
phi2 = phi21(k)+ offset2;

end