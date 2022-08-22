function sol = FreeRK4SE3N(f,action,h,p)
        
    %Updates the solution p one time step in the future with RKMK method of
    %order 4 built commutator free.

    gAc = @(g,x) action(exponentialSE3N(g),x); 
    
    Y1 = p;
    k1 = h*f(Y1);
    
    Y2 = gAc(k1/2,p);
    k2 = h*f(Y2);
    
    Y3 = gAc(k2/2,p);
    k3 = h*f(Y3);
    
    Y4 = gAc(k3-k1/2,Y2);
    k4 = h*f(Y4);
    
    yHalf = gAc(1/12*(3*k1+2*k2+2*k3-k4),p);
    sol = gAc(1/12*(-k1+2*k2+2*k3+3*k4),yHalf); 
end