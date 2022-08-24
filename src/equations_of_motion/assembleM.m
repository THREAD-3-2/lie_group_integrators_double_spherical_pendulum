function M = assembleM(q, L, m)
%
% :param q: position vector [q1, ..., qN] in S^2
% :param L: length of the pendulum
% :param m: mass of the pendulum
%
% :returns: inertia mass matrix of the system

    N = floor(length(q) / 3); % Number of connected pendulums
    
    M = zeros(3 * N, 3 * N);

    for i = 1 : N
        for j = i + 1 : N
            M(3 * i - 2 : 3 * i, 3 * j - 2 : 3 * j) = sum(m(j : end)) * L(i) * L(j) * eye(3);
        end
    end
    
    M = M + M';
    
    for i = 1 : N
        M(3 * i - 2 : 3 * i, 3 * i - 2 : 3 * i) = sum(m(i : end)) * L(i) ^ 2 * eye(3);
    end
end
