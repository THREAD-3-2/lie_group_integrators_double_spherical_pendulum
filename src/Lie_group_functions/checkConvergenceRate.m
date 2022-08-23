function [] = checkConvergenceRate(f,action,vecField,z0,L,m)

    T = 0.1;
    
    % To highlight the right convergence rate, it is important to select
    % accurately the number of Timestep to consider, and this may vary
    % accordingly to the initial condition.
    
    n = 2 : 5;
    N0 = 5;
    Nrange = N0 * 2.^n;
    
    %Define the vectors where we store the errors wrt a reference solution with ODE45
    errorLE1 = zeros(length(Nrange),1);
    errorEH1 = errorLE1;
    errorRK4FR1 = errorLE1;
    errorRK4_2Comm1 = errorLE1;
    errorRK41 = errorLE1;
    errorRK31 = errorLE1;

    %Solve with ODE45
    q0 = extractq(z0);
    w0 = extractw(z0);
    
    z0Ref = [q0; w0]; %We reorder the initial condition to make ode45 work withour reordering the equations i.e. we pass from [q_1,w_1,...,q_N,w_N] to [q_1,q_2,...,q_N,w_1,w_2,...,w_N]
    odeFunc = @(t, z) [FuncQ(z); FuncW(z, L, m)]; %FuncQ,FuncW just assemble the first N and last N equations
    options = odeset('AbsTol', 1e - 12, 'RelTol', 1e - 12);
    [timeSol,zSol] = ode45(odeFunc, [0 T], z0Ref, options);
    zSol = zSol';
    zSol = reorder(zSol); %we pass from [q_1,q_2,...,q_N,w_1,w_2,...,w_N] to [q_1,w_1,...,q_N,w_N] so that we can compare it easily with the
    %solution coming from our numerical schemes
    zRef = zSol(:, end);
    
    %We compute even a reference solution with Commutator free method
    N = 2000;
    zRef4 = z0;
    time = linspace(0, T, N);
    dt = time(2) - time(1);
        
    for i = 1 : N - 1
        zRef4 = COMMFREERK4(f, action, dt, zRef4);
    end

    count = 1;
    dts = zeros(length(Nrange), 1);
    
    %Solve in time for the different Timestep and compute the errors
    for N = Nrange
        
        time = linspace(0, T, N);
        dt = time(2) - time(1);
        dts(count) = dt;
        z1 = z0;
        z2 = z0;
        z3 = z0;
        z4 = z0;
        z5 = z0;
        z6 = z0;
        
        for i = 1 : N - 1
            z1 = LieEulerSE3N(vecField, action, z1, dt);
            z2 = LieEulerHeunSE3N(vecField, action, z2, dt);
            z3 = COMMFREERK4(f, action, dt, z3);
            z4 = RKMK4(vecField, action, z4, dt);
            z5 = TwoCommRK4SE3N(f, action, dt, z5);
            z6 = LieRK3_SE3N(vecField, action, z6, dt);
        end

        %% Comparison with ODE45
        errorLE1(count) = norm(z1 - zRef, 2);
        errorEH1(count) = norm(z2 - zRef, 2);
        errorRK4FR1(count) = norm(z3 - zRef, 2);
        errorRK41(count) = norm(z4 - zRef, 2);
        errorRK4_2Comm1(count) = norm(z5 - zRef, 2);
        errorRK31(count) = norm(z6 - zRef, 2);
        
        count = count + 1;
    end

    figure;
    %% Convergence rate of Lie Euler
    ord1 = (dts/dts(1)).^(1) * errorLE1(1); 
    ord2 = (dts/dts(1)).^(2) * errorLE1(1);
    ord3 = (dts/dts(1)).^(3) * errorLE1(1);
    ord4 = (dts/dts(1)).^(4) * errorLE1(1);

    subplot(2, 3, 1) 
    loglog(dts, errorLE1, 'k-*', dts, ord1, 'r--', dts, ord2, 'c--', dts, ord3, 'm--', dts, ord4, 'g--', 'LineWidth', 2)
    xlabel('Timestep')
    ylabel('Global error')
    title("Lie Euler",'FontSize', 14)
    legend('Global error', 'ord1', 'ord2', 'ord3', 'ord4', 'FontSize', 14, 'Location', 'Best');
    ax = gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;

    %% Convergence rate of Euler Heun
    ord1 = (dts/dts(1)).^(1) * errorEH1(1); 
    ord2 = (dts/dts(1)).^(2) * errorEH1(1);
    ord3 = (dts/dts(1)).^(3) * errorEH1(1);
    ord4 = (dts/dts(1)).^(4) * errorEH1(1);

    subplot(2, 3, 2)
    loglog(dts, errorEH1, 'k-*', dts, ord1, 'r--', dts, ord2, 'c--', dts, ord3, 'm--', dts, ord4, 'g--', 'LineWidth', 2)
    xlabel('Timestep')
    ylabel('Global error')
    title("Euler Heun",'FontSize', 14)
    legend('Global error', 'ord1', 'ord2', 'ord3', 'ord4', 'FontSize', 14, 'Location', 'Best');
    ax = gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;


    %% Convergence rate of RKMK4
    ord1 = (dts/dts(1)).^(1) * errorRK41(1);
    ord2 = (dts/dts(1)).^(2) * errorRK41(1);
    ord3 = (dts/dts(1)).^(3) * errorRK41(1);
    ord4 = (dts/dts(1)).^(4) * errorRK41(1);

    subplot(2,3, 3)
    loglog(dts, errorRK41, 'k-*', dts, ord1, 'r--', dts, ord2, 'c--', dts, ord3, 'm--', dts, ord4, 'g--', 'LineWidth', 2)
    xlabel('Timestep')
    ylabel('Global error')
    title("RKMK4",'FontSize', 14)
    legend('Global error', 'ord1', 'ord2', 'ord3', 'ord4', 'FontSize', 14, 'Location', 'Best');
    ax = gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;

     %% Convergence rate of RKMK4 commutator Free
    ord1 = (dts/dts(1)).^(1) * errorRK4FR1(1);
    ord2 = (dts/dts(1)).^(2) * errorRK4FR1(1);
    ord3 = (dts/dts(1)).^(3) * errorRK4FR1(1);
    ord4 = (dts/dts(1)).^(4) * errorRK4FR1(1);

    subplot(2,3, 4)
    loglog(dts, errorRK4FR1, 'k-*', dts, ord1, 'r--', dts, ord2, 'c--', dts, ord3, 'm--', dts, ord4, 'g--', 'LineWidth', 2)
    xlabel('Timestep')
    ylabel('Global error')
    title("RKMK4 Commutator FREE",'FontSize', 14)
    legend('Global error', 'ord1', 'ord2', 'ord3', 'ord4', 'FontSize', 14, 'Location', 'Best');
    ax = gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;

    %% Convergence rate of RKMK4 2 commutators
    ord1 = (dts/dts(1)).^(1) * errorRK4_2Comm1(1);
    ord2 = (dts/dts(1)).^(2) * errorRK4_2Comm1(1);
    ord3 = (dts/dts(1)).^(3) * errorRK4_2Comm1(1);
    ord4 = (dts/dts(1)).^(4) * errorRK4_2Comm1(1);

    subplot(2, 3, 5)
    loglog(dts, errorRK4_2Comm1, 'k-*', dts, ord1, 'r--', dts, ord2, 'c--', dts, ord3, 'm--', dts, ord4, 'g--', 'LineWidth', 2)
    xlabel('Timestep')
    ylabel('Global error')
    title("RKMK4 with 2 Commutators",'FontSize', 14)
    legend('Global error', 'ord1', 'ord2', 'ord3', 'ord4', 'FontSize', 14, 'Location', 'Best');
    ax = gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;

    %% Convergence rate of RKMK3 
    ord1 = (dts/dts(1)).^(1) * errorRK31(1);
    ord2 = (dts/dts(1)).^(2) * errorRK31(1);
    ord3 = (dts/dts(1)).^(3) * errorRK31(1);
    ord4 = (dts/dts(1)).^(4) * errorRK31(1);

    subplot(2, 3, 6)
    loglog(dts, errorRK31, 'k-*', dts, ord1, 'r--', dts, ord2, 'c--', dts, ord3, 'm--', dts, ord4, 'g--', 'LineWidth', 2)
    xlabel('Timestep')
    ylabel('Global error')
    title("RKMK3",'FontSize', 14)
    legend('Global error', 'ord1', 'ord2', 'ord3', 'ord4', 'FontSize', 14, 'Location', 'Best');
    ax = gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;

    sgtitle("Convergence rate comparison against ODE45",'FontSize', 16);
    set(gcf, 'Position',  [0   200   500   700])
    
end
