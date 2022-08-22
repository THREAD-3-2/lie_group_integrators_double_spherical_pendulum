function sol =  RKMK4(vector_field,exponential,action,p,h,sigma0,trajectory,t)
% Runge-Kutta-Munthe-Kaas time integrator order 4
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
% :returns: solution at time t_(n+1)

    k1 = vector_field(sigma0,p,trajectory(t));
    k2 = vector_field(h/2*k1, p,trajectory(t+h/2));
    k3 = vector_field(h/2 * k2 , p,trajectory(t+h/2));
    k4 = vector_field(h * k3 , p,trajectory(t+h));

    sol = action(exponential(h/6*k1 + h/3*k2 + h/3*k3 + h/6*k4),p);
       
end
