function A = expSE3(input)
% Exponential map on SE(3)
%
% :param input: element of the Lie algebra se(3) represented as [hat(u),v],
%               with u and v 3x1 column vectors
%
% :returns: element of the group SE(3)
    
    u = input(1:3);
    v = input(4:6);
    theta = norm(u,2);
    
    tol = 1e - 16;  
    
    if theta > tol
        A = sin(theta) / theta;
        B = (1 - cos(theta)) / (theta^2);
        C = (1 - A)/(theta ^ 2);
        V = eye(3) + B * hat(u) + C * hat(u) * hat(u);
        A = [expRodrigues(u), V * v];
    elseif theta == 0
        A = [expRodrigues(u), v];
    else
        Blow = 0.5 - theta^2/24 + theta^4/720 - theta^6/40320;
        Clow = (1/6 - theta^2/120 + theta^4/5040 - theta^6/362880);
        Vlow = eye(3) + Blow * hat(u) + Clow * hat(u) * hat(u);
        A = [expRodrigues(u), Vlow * v];
    end
end
