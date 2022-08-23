function [v, N] = tangentBehaviour(f, action, vecField, z0, L, m)


    %In this method we solve the equation and compute the tangentiality
    %conditions for the solutions obtained with the various schemes. The
    %output value N is just the number of time instants considered.

    T = 5; %Final time
    
    M = length(z0)/6; %Number of connected pendulums

    q0 = extractq(z0);
    w0 = extractw(z0);
    
    z0Ref = [q0; w0];
    N = 1000;
    time = linspace(0, T, N);
    dt = time(2) - time(1);
    odeFunc = @(t, z) [FuncQ(z); FuncW(z, L, m)]; %FuncQ,FuncW just assemble the first N and last N equations
    options = odeset('AbsTol', 1e - 12, 'RelTol', 1e - 12);
    [timeSol,zSol] = ode45(odeFunc, time, z0Ref, options); %Solve it with ODE45
    zSolODE45 = zSol';
    zSolODE45 = reorder(zSolODE45);

    zSol = RK4(odeFunc, time, z0Ref); %Solve it with Runge Kutta 4
    zSol = reorder(zSol);
    
    v = zeros(M, N, 8); %Tensor where we store the norms
    v(:, :, 1) = checkTangency(zSolODE45);
    v(:, :, 2) = checkTangency(zSol);
            
    z1 = zeros(6 * M, N);
    z2 = z1;
    z3 = z1;
    z4 = z1;
    z5 = z1;
    z6 = z1;
    z1(:, 1) = z0;
    z2(:, 1) = z0;
    z3(:, 1) = z0;
    z4(:, 1) = z0;
    z5(:, 1) = z0;
    z6(:, 1) = z0;

    for i = 1 : N - 1
        z1(:, i + 1) = LieEuler(vecField, action, z1(:, i), dt);
        z2(:, i + 1) = LieEulerHeun(vecField, action, z2(:, i), dt);
        z3(:, i + 1) = COMMFREERK4(f, action, dt, z3(:, i));
        z4(:, i + 1) = RKMK4(vecField, action, z4(:, i), dt);
        z5(:, i + 1) = TwoCommRKMK4(f, action, dt, z5(:, i));
        z6(:, i + 1) = RKMK3(vecField, action, z6(:, i), dt);
    end
    
    v(:, :, 3) = checkTangency(z1);
    v(:, :, 4) = checkTangency(z2);
    v(:, :, 5) = checkTangency(z3);
    v(:, :, 6) = checkTangency(z4);
    v(:, :, 7) = checkTangency(z5);
    v(:, :, 8) = checkTangency(z6);

end