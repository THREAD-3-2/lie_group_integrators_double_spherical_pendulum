function sol =  LieRK4_SE3N(vecField,action,p,h)
    

    %Computes one time step update with RKMK of order 4

    sigma0 = zeros(length(p),1);
    k1 = vecField(sigma0,p);
    k2 = vecField(h/2*k1, p);
    k3 = vecField(h/2 * k2 , p);
    k4 = vecField(h * k3 , p);


    sol = action(exponentialSE3N(h/6*k1 + h/3*k2 + h/3*k3 + h/6*k4),p);      
    
       
end
