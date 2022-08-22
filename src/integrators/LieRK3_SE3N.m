function sol =  LieRK3_SE3N(vecField,action,p,h)
    
    %Computes one time step update with RKMK method of order 3

    k1 = vecField(zeros(length(p),1),p);
    k2 = vecField(h/2*k1,p);
    k3 = vecField(h*(-k1+2*k2), p);
    sol = action(exponentialSE3N(h*(k1/6+2*k2/3+k3/6)),p);
       
end