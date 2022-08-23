function sol = LieEuler(vecField, action, p, h)
% Lie-Euler time integrator (Runge-Kutta-Munthe-Kaas method of order 1)
%
% :param vecField: right hand side of the ODE
% :param action: Lie group action
% :param p: solution at time t_n
% :param h: time step size
%
% :returns: solution at time t_(n+1) order 1

    k0 = zeros(length(p), 1);
    k1 = vecField(k0, p);

    sol = action(exponentialSE3N(h * k1), p);   
    
end
