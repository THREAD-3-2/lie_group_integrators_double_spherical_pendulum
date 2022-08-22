function res = exponentialSE3N(sigma)
% Exponential map on SE(3)^N
%
% :param sigma: element of se(3)^N
%
% :returns: element of SE(3)^N


    N = floor(length(sigma)/6);
    res = zeros(3 * N, 4);
    for i = 1 : N
       res(3 * i - 2 : 3 * i, :) = expSE3(sigma(6 * i - 5 : 6 * i));
    end
end
