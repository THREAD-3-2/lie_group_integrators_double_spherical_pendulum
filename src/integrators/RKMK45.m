function [sol, err] = variableRKMK45(vecField, exponentialSE3N, action, p, h, sigma0)
% Runge-Kutta-Munthe-Kaas pair (5,4) based on the Dormand-Prince method DOPRI5(4)
%
% :param vecField: right hand side of the ODE
% :param exponentialSE3N: exponential map from se(3)^N to SE(3)^N
% :param action: Lie group action
% :param p: solution at time t_n
% :param h: time step size
% :param sigma0: initial value of the curve sigma on se(3)^N
%
% :returns: solution at time t_(n+1) and local error

    % The one step calculation in the Dormand-Prince method DOPRI5(4) is done as follows:
    sigma0 = zeros(length(p), 1);
    k1 = vecField(sigma0, p);
    k2 = vecField(h / 5 * k1, p);
    k3 = vecField(h * (3 / 40 * k1 + 9 / 40 * k2), p);
    k4 = vecField(h * (44 / 45 * k1 - 56 / 15 * k2 + 32 / 9 * k3), p);
    k5 = vecField(h * (19372 / 6561 * k1 - 25360 / 2187 * k2 + 64448 / 6561 * k3 - 212 / 729 * k4), p);
    k6 = vecField(h * (9017 / 3168 * k1 - 355 / 33 * k2 + 46732 / 5247 * k3 + 49 / 176 * k4 - 5103 / 18656 * k5), p);
    k7 = vecField(h * (35 / 384 * k1 + 500 / 1113 * k3 + 125 / 192 * k4 - 2187 /6784 * k5 + 11 / 84 * k6), p);
    
    % Then the next step value y_{n+1} is calculated by Runge-Kutta method of order 4
    sigma = h * (35 / 384 * k1 + 500 / 1113 * k3 + 125 / 192 * k4 - 2187 / 6784 * k5 + 11 / 84 * k6);
    
    % Next, we will calculate the next step value z_{k+1} by Runge-Kutta method of order 5 
    sigmaHat = h * (5179 / 57600 * k1 + 7571 / 16695 * k3 + 393 / 640 * k4 - 92097 / 339200 * k5 + 187 / 2100 * k6 + 1 / 40 * k7);
            
    sol = action(exponentialSE3N(sigma), p);  
            
    % We take the difference at the level of the Lie algebra approximations
    err = norm(sigma - sigmaHat, 2); %We take the difference at the level of the Lie algebra approximations
       
end
