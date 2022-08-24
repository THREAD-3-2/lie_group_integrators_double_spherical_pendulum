function vec = FuncQ(z)
% This function is used to integrate with ODE45
%
% :param z: is of the form z = [q1, q2, ..., qP, w1, w2, ..., wP]
%
% :returns: the part of the vector field for the \dot{q}_i, so hat(w_i)*q_i

    q = z(1 : length(z)/2);
    w = z(length(z)/2 + 1 : end);
    
    vec = zeros(length(z)/2, 1);
    P = length(z)/6;
    
    for i = 1 : P
        vec(3 * i - 2 : 3 * i) = cross(w(3 * i - 2 : 3 * i), q(3 * i - 2 : 3 * i));
    end
end
