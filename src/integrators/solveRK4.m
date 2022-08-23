function z = RK4(f, time, z0)
% Classical implementation of Runge-Kutta of order 4
%
% :param f: right hand side of the ODE
% :param time: time
% :param z0: initial value
%
% :returns: solution at time t_(n+1)

    N = length(time); 
    z = zeros(length(z0), N);
    z(:, 1) = z0;
    
    for i = 1 : N - 1
        
        dt = time(i + 1) - time(i);
        t = time(i);
        p = z(:, i);
        
        k1 = dt * f(t, p);
        k2 = dt * f(t + dt/2, p + k1/2);
        k3 = dt * f(t + dt/2, p + k2/2);
        k4 = dt * f(t + dt, p + k3);
        
        z(:, i + 1) = z(:, i) + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4); 
    end
end
