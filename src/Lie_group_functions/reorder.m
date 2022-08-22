function vec = reorder(z)

    %This method turns the vector z = [q1,q2,..,qN,w1,w2,...,wN] into
    % vec = [q1,w1,q2,w2,...,qN,wN] which is the format used in the Lie
    % group integrators setting.

    N = length(z(:, 1))/6; %Number of connected pendulums
    q = z(1 : length(z(:, 1))/2, :);
    w = z(length(z(:, 1))/2 + 1 : end, :);
    
    vec = zeros(size(z));
    
    for i = 1 : 2 * N
        if mod(i, 2) == 0 
            vec(3 * i - 2 : 3 * i, :) = w(3 * i/2 - 2 : 3 * i/2, :);
        else
            vec(3 * i - 2 : 3 * i, :) = q(3 * (i + 1)/2 - 2 : 3 * (i + 1)/2, :);
        end
    end
end