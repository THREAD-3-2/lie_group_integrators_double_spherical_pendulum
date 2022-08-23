function sol = CommFreeRK4(f, action, h, p)
% Commutator-free time integrator of order 4
%
% :param f: map f from the phase space (on which the vector field is defined) to the Lie algebra
% :param action: Lie group action
% :param h: time step size
% :param p: solution at time t_n
%
% :returns: solution at time t_(n + 1)

    gAc = @(g,x) action(exponentialSE3N(g),x); 
    
    Y1 = p;
    k1 = h * f(Y1);
    
    Y2 = gAc(k1/2, p);
    k2 = h * f(Y2);
    
    Y3 = gAc(k2/2, p);
    k3 = h * f(Y3);
    
    Y4 = gAc(k3 - k1/2, Y2);
    k4 = h * f(Y4);
    
    yHalf = gAc(1/12 * (3 * k1 + 2 * k2 + 2 * k3 - k4), p);
    sol = gAc(1/12 * (- k1 + 2 * k2 + 2 * k3 + 3 * k4), yHalf); 
end
