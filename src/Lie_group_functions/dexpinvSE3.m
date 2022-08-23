function vec = dexpinvSE3 (sigma, input)
% Inverse of the derivative of the exponential map on SE(3)
%
% :param sigma: element in the Lie algebra se(3), represented as a 6x1 vector  
% :param input: element in the Lie algebra se(3), represented as a 6x1 vector  
%
% :returns: inverse of dexp_sigma(input) as 6x1 vector
    
    A = sigma(1:3);
    a = sigma(4:6);
    alpha = norm(A, 2);
    
    rho = A' * a;
    g1 =  - 0.5;
    g1tilde =  0;
    f0 = 1;
    mat = zeros(3);
    
    tol = 1e-16;

    if (alpha) > tol
        g2 = (1 - alpha/2 * cot(alpha/2))/(alpha^2);
        g2prime = 1/alpha^2 * (- 0.5 * cot(alpha/2) - alpha/2 * (- 0.5 - 0.5 * cot(alpha/2)^2)) - 2/alpha^3 * (1 - alpha/2 * cot(alpha/2));
        g2tilde = rho/alpha * g2prime;
        
        Mat = f0 * eye(6) + [g1 * hat(A), mat; g1 * hat(a) + g1tilde * hat(A), g1 * hat(A)] + [hat(A), mat; hat(a), hat(A)] * [g2 * hat(A), mat; g2 * hat(a) + g2tilde * hat(A), g2 * hat(A)];
        vec = Mat * input;
    else
        g2low = 1/12 + alpha^2/720 + alpha^4/(30240) + alpha^6/(1209600);
        g2tildelow = rho * (1/360 + alpha^2/7560);
        
        Mat = f0 * eye(6) + [g1 * hat(A), mat; g1 * hat(a) + g1tilde * hat(A), g1 * hat(A)] + [hat(A), mat; hat(a), hat(A)] * [g2low * hat(A), mat; g2low * hat(a) + g2tildelow * hat(A), g2low * hat(A)];
        vec = Mat * input;
    end
end
