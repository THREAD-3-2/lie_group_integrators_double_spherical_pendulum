function [q0, w0, z0] = initializeStat(N)
% Initialization of the system
%
% :param N: number of connected pendulums
%
% :returns: initial positions for q and w of the N-fold pendulum

    z0 = zeros(6 * N, 1);
    qref = [1; 0; 0];
    w0 = 0 * qref;
    
    for i = 1 : N
        z0(6 * i - 5 : 6 * i) = [qref; w0];
    end
    
    q0 = extractq(z0);
    w0 = extractw(z0);
        
end
