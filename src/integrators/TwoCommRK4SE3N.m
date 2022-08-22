function sol = TwoCommRK4SE3N(f,action,h,p)
        
    %RKMK 4 with the number of commutators reduced to 2

    gAc = @(g,x) action(exponentialSE3N(g),x); 
    
    Y1 = p;
    k1 = h*f(Y1);
    
    Y2 = gAc(k1/2,p);
    k2 = h*f(Y2);
    
    Y3 = gAc(k2/2- 1/8*commutatorSE3N(k1,k2),p);
    k3 = h*f(Y3);
    
    Y4 = gAc(k3,p);
    k4 = h*f(Y4);
    
    sol = gAc(1/6*(k1+2*k2+2*k3+k4-0.5*commutatorSE3N(k1,k4)),p); 

end