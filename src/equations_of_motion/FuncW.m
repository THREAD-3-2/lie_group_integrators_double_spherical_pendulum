function vec = FuncW(z,L,m)
% This function is used to integrate with ODE45
%
% :param z: is of the form z = [q1, q2, ..., qN, w1, w2, ..., wN]
% :param L: length of the pendulum
% :param m: mass of the pendulum
%
% :returns: the part of the vector field for the \dot{w}_i, so R(q)^{-1}*RHS of the ODE, defined by assembleF

    q = z(1 : length(z)/2);
    w = z(length(z)/2 + 1 : end);
    
    R = assembleR(q, L, m);
    F = assembleF(q, w, m, L);
    
    vec = R\F;
end
