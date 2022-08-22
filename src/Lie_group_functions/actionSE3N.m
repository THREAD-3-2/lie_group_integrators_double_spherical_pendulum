function res = actionSE3N(B,input)
% action on SE(3)^N
%
% :param B: element of SE(3)^N
% :param input: element of se(3)^N
%
% :returns: element of se(3)^N

    N = length(input)/6;   % number of connected pendulums
    
    res = zeros(6 * N, 1);
    
    for i = 1:N
        res(6 * i - 5 : 6 * i) = actionSE3(B(3 * i - 2 : 3 * i, :), input(6 * i - 5 : 6 * i));  % Assemble the various actions of SE(3)
    end
end
