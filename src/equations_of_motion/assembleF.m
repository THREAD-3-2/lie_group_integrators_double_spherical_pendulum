function F = assembleF(q, w, m, L)
% This function defines the right hand side, precisely we need it to define the equations for the angular velocities, which now becomes R(q)w' = F, and here we assemble this F vector.
%
% :param q: position vector [q1, ..., qN] in S^2
% :param w: angular velocities vector [w1, ..., wN] in T_{qi}S^2
% :param L: length of the pendulum
% :param m: mass of the pendulum
%
% :returns: equations for the angular velocities of the form R(q)w' = F

    N = length(m); % Number of connected pendulums
    f = cell(1, N);
    g = 9.81;
    e3 = [0; 0; 1];
    
    M = assembleM(q, L, m);
    
    Block = @(Mat, i, j) getBlock(Mat, i, j);
    vec = @(v, i) getVec(v, i);

    for i = 1 : N
        f{i} = - hat(vec(q, i)) * sum(m(i : end)) * g * L(i) * e3;
        for j = 1 : N
            if j ~= i
                f{i} = f{i} + Block(M, i, j) * norm(vec(w, j), 2) ^ 2 * hat(vec(q, i)) * vec(q, j);
            end
        end
    end
    
    F = zeros(3 * N, 1);
    for i = 1 : N
        F(3 * i - 2 : 3 * i) = f{i};
    end

end
