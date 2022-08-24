function q = extractq(v)
% extracts the q components from [q, w]
%
% :param v: input vector v = [q1, w1, q2, w2, ..., qN, wN]
%
% :returns: the q components of v, q = [q1,q2,...,qN]

    N = length(v(:, 1))/6;   % Number of connected pendulums
    m = length(v(1, :));
    q = zeros(3 * N, m);
    for i = 1 : N
        q(3 * i - 2 : 3 * i, :) = v(6 * i - 5 : 6 * i - 3, :);
    end
end
