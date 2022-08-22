function w = extractw(v)
% extracts the w components from [q, w]
%
% :param v: input vector v = [q1, w1, q2, w2, ..., qN, wN]
%
% :returns: the w components of v, w = [w1, w2, ..., wN]

    N = length(v)/6;  %Number of connected pendulums
    w = zeros(3 * N, 1);
    for i = 1 : N
        w(3 * i - 2 : 3 * i) = v(6 * i - 2 : 6 * i);
    end
end
