function res = dexpinvSE3N(sigma, input)
% Inverse of the derivative of the exponential map on SE(3)^N
%
% :param sigma: element of se(3)^N 
% :param input: element of se(3)^N 
%
% :returns: inverse of dexp_sigma(input)

    N = floor(length(sigma)/6);
    res = zeros(6 * N, 1);
    for i = 1 : N
       res(6 * i - 5 : 6 * i) = dexpinvSE3(sigma(6 * i - 5 : 6 * i), input(6 * i - 5 : 6 * i));
    end
end
