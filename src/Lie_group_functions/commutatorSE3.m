function res = commutatorSE3(u,v)
    
    %This method computes the commutator [u,v] = [(u1,u2),(v1,v2)] = 
    %(u1 x v1, u1 x v2 - v1 x u2).

    u1 = u(1:3);
    u2 = u(4:end);
    
    v1 = v(1:3);
    v2 = v(4:end);
    
    res = zeros(6,1);
    res(1:3) = hat(u1)*v1;
    res(4:end) = hat(u1)*v2 - hat(v1)*u2;
    

end