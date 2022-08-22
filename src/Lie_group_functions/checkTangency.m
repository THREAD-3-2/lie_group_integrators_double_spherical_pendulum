function v = checkTangency(z)

    %As an input z we get the matrix with the solutions, each column is a
    %set [q1,w1,...,qN,wN] at a specific time t_j, where j is the column of
    %v where this solution is stored.

    N = length(z(:, 1))/6; %Number of connected pendulums
    l = length(z(1, :)); %Number of time steps
    
    v = zeros(N, l); %Each colum is associated to a time instant, each row to a specific pendulum
    %v(i, j) will be the the tangentiality contidion at time tj for the i-th
    %pendulum 
    
    for i = 1 : N
        for j = 1 : l
            v(i, j) = z(6 * i - 5 : 6 * i - 3, j)' * z(6 * i - 2 : 6 * i, j);
        end
    end
end