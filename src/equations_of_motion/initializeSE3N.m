function [q0, w0, z0] = initializeSE3N(N)
% This method randomly picks a point in (TS^2)^N
%
% :param N: number of the connected pendula
%
% :returns: randomly picked point in (TS^2)^N

    w0 = rand(3 * N, 1);
    q0 = w0;
    z0 = rand(6 * N, 1);
    
    for i = 1 : N
        q0(3 * i - 2 : 3 * i) = rand(3, 1);
        q0(3 * i - 2 : 3 * i) = q0(3 * i - 2 : 3 * i)/norm(q0(3 * i - 2 : 3 * i), 2);
        v = rand(3, 1);
        w0(3 * i - 2 : 3 * i) = cross(q0(3 * i - 2 : 3 * i), v);
        z0(6 * i - 5 : 6 * i) = [q0(3 * i - 2 : 3 * i); w0(3 * i - 2 : 3 * i)];
    end
end
