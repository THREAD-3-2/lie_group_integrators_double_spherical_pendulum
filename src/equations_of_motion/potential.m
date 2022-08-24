function val = potential(q, L, m)
%
% :param q: position vector [q1, ..., qN] in S^2
% :param L: length of the pendulum
% :param m: mass of the pendulum
%
% :returns: potential energy of the system

    P = length(L);
    g = 9.81;
    e3 = [0; 0; 1];
    
    val = 0;
    for i = 1 : P
        val = val + L(i) * g * sum(m(i : end)) * e3' * q(3 * i - 2 : 3 * i);
    end
end
