function res = commutatorSE3N(u, v)
% Computes the commutator
%
% :param u: element of se(3)^N 
% :param v: element of se(3)^N 
%
% :returns: commutator [u, v]
    
    N =  length(u)/6;    %number of connected pendulums
    res = zeros(length(u), 1);
    
    for i = 1 : N
        res(6 * i - 5 : 6 * i) = commutatorSE3(u(6 * i - 5 : 6 * i), v(6 * i - 5 : 6 * i)); %exploit the direct product structure of the Lie group SE(3)^N
    end
end
