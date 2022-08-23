function sol = LieEulerHeunSE3N(vecField, action, p, h)
% Runge-Kutta-Munthe-Kaas time integrator order 2, based on Heun scheme
%
% :param vecField: right hand side of the ODE
% :param action: Lie group action
% :param p: solution at time t_n
% :param h: time step size
%
% :returns: solution at time t_(n+1)

    k0 = zeros(length(p), 1);
    k1 = vecField(k0, p);
    k2 = vecField(h * k1, p);
    sol = action(exponentialSE3N(h/2 * (k1 + k2)), p);

end
