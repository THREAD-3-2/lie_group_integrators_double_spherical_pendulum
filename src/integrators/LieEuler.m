function sol = LieEuler(vector_field,exponential,action,p,h,sigma0,trajectory,t)
% Lie-Euler time integrator (Runge-Kutta-Munthe-Kaas method of order 1)
%
% :param vector_field: right hand side of the ODE
% :param action: Lie group action
% :param expinential: exponential map from the Lie algebra to the Lie group
% :param p: solution at time t_n
% :param h: time step size
% :param sigma0: initial value of the curve sigma on the Lie algebra
% :param trajectory: desired trajectory
% :param t: discrete time t_n
%
% :returns: solution at time t_(n+1))order 1

    sol = action(exponential(h*vector_field(sigma0,p,trajectory(t))),p);
end