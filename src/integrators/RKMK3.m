function sol =  RKMK3(vecField, action, p, h)
% Runge-Kutta-Munthe-Kaas time integrator order 3
%
% :param vecField: right hand side of the ODE
% :param action: Lie group action
% :param p: solution at time t_n
% :param h: time step size
%
% :returns: solution at time t_(n + 1)

    k1 = vecField(zeros(length(p), 1), p);
    k2 = vecField(h/2 * k1, p);
    k3 = vecField(h * (- k1 + 2 * k2), p);
    sol = action(exponentialSE3N(h * (k1/6 + 2 * k2/3 + k3/6)), p);
       
end
