function sol =  RKMK4(vecField, action, p, h)
% Runge-Kutta-Munthe-Kaas time integrator order 4
%
% :param vecField: right hand side of the ODE
% :param action: Lie group action
% :param p: solution at time t_n
% :param h: time step size
%
% :returns: solution at time t_(n+1)

    sigma0 = zeros(length(p), 1);
    k1 = vecField(sigma0, p);
    k2 = vecField(h/2 * k1, p);
    k3 = vecField(h/2 * k2 , p);
    k4 = vecField(h * k3 , p);

    sol = action(exponentialSE3N(h/6 * k1 + h/3 * k2 + h/3 * k3 + h/6 * k4), p);      
       
end
