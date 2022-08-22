% main function
%
%This code contains the numerical experiments with the step adaptive RKMK method based on Dormand-Prince pair (5,4)
%
% Ref.: E. Celledoni, E. Çokaj, A. Leone, D. Murari, B. Owren.
% "Dynamics of the N-fold Pendulum in the framework of Lie Group Integrators", arXiv.
%
% Ref.: E. Celledoni, E. Çokaj, A. Leone, D. Murari, B. Owren.
% "Lie group integrators for mechanical systems",
% International Journal of Computer Mathematics, 99:1, 58-88.


clc;
clear all;
close all;

%% Setting the parameters 

Prange = 2 : 2 : 20;    %number of connected pendulums
accVar = zeros(length(Prange), 1);
accConst = accVar;
index = 1;
tol = 1e - 6;
Lref = 5;               % length of the entire chain of pendulum
steps = zeros(length(Prange), 1);

for P = Prange
    
    L = rand(P, 1) + 0.5; % Set random lengths
    m = rand(P, 1) + 0.5; % Set random masses
    L = 0 * L + Lref / P; % size of each pendulum
    m = 0 * m + 1;        % mass of the pendulum
 
    [q0, w0, z0] = initializeStat(P);  % initial positions and angular velocities of the N-fold pendulum

    t0 = 0;         % initial time
    T = 3;          % final time
    N = 1000;       % number of time steps
    time = linspace(t0, T, N); 
    dt = time(2) - time(1);   % stepsize

    getq = @(v) extractq(v);
    getw = @(v) extractw(v);

    f = @(v) fManiToAlgebra(getq(v), getw(v), L, m); 
    action = @(B, input) actionSE3N(B, input); 
    vecField = @(sigma, p) dexpinvSE3N(sigma, f(action(exponentialSE3N(sigma), p)));

    z = z0;
    qC = zeros(3 * P,N);
    pC = qC;

    Len = zeros(3 * P, 1);
    for i = 1 : P
        Len(3 * i - 2 : 3 * i) = L(i) * ones(3, 1);
    end
 
    Mat = diag(Len);
    if P > 1
        for i = 3 : 3 : 3 * (P - 1)
            Mat = Mat + diag(Len(1 : 3 * P - i), -i);
        end
    end
    qC(:, 1) = q0;
    pC(:, 1) = Mat * q0;


    %% COMPARISON BETWEEN RKMK45 AND RKMK5, i.e. VARIABLE STEPSIZE AGAINST CONSTANT ONE

    chunk = 400;
    Z = zeros(length(z0), chunk);
    TT = zeros(1, chunk);
    Y = Z;
    Y(:, 1) = z0;

    a = 1/4;
    theta = 0.85;
    i = 1;
    h = T / 100;
    rejected = 0;
    ctr = 0;
    
    while TT(i) < T - 5 * eps
        err = tol + 1;
        while err > tol 
            [z, err] = variableRKMK45(vecField, action, Y(:, i), h);
            accepted = (err < tol);
            if accepted
                i = i + 1;
                Y(:, i) = z;
                TT(i) = TT(i - 1) + h;

                if ctr == 0
                    ctr = 1;
                    varStrict = z;
                end            
            else
                rejected = rejected + 1;
            end
            h = min(theta * (tol / err) ^ a * h, T - TT(i));
        end
        if mod(i, chunk) == 0
            Y = [Y, Z];
            TT = [TT zeros(1, chunk)];
        end
    end
    Y = Y(:, 1 : i);
    TT = TT(:, 1 : i);

    Nsteps = i - 1;
    tt = linspace(0, T, Nsteps + 1);
    dt = tt(2) - tt(1);
    zRK45 = zeros(length(z0), Nsteps);
    zRK45(:, 1) = z0;
    count = 1;
    ts = 0;

    for k = 1 : Nsteps
        zRK45(:, k + 1) = RKMK5(vecField, action, zRK45(:, k), dt);
        ts = ts + dt;
    end
     
    % Solve with ODE45

    z0Ref = [q0; w0];
    odeFunc = @(t, z) [FuncQ(z); FuncW(z, L, m)];
    options = odeset('AbsTol', 1e - 12, 'RelTol',1e - 6);
    [timeSol, zC] = ode45(odeFunc, [0 T], z0Ref, options);
    zC = zC';
    zC = reorder(zC);
    zRef = zC(:, end);

    accVar(index) = norm(zRef - Y(:, end), 2);
    accConst(index) = norm(zRef - zRK45(:, end), 2);
    steps(index) = Nsteps;
    index = index + 1;
end

figure;
yyaxis left
semilogy(Prange, accConst, 'r-o', 'linewidth', 3, 'Markersize', 10);
hold on;
semilogy(Prange, accVar, 'k-*', 'linewidth', 3, 'Markersize', 10);
xlabel('Number of connected pendulums', 'FontSize', 30);
ylabel('Accuracy at T = 3', 'FontSize', 30);
hold on;
ax = gca;
set(ax, 'xcolor', 'k')
set(ax, 'ycolor', 'k')
yyaxis right
ax = gca;
plot(Prange, steps, 'b-s', 'linewidth', 3, 'Markersize', 10);
set(ax, 'ycolor', 'b')
ax.FontSize = 20;
ylabel('Number of time steps', 'FontSize', 30);
legend('Accuracy RKMK5', 'Accuracy RKMK(5, 4)', 'Number of time steps', 'FontSize', 30);

%% Comparison with step adaptation in ODE45
figure;
z0Ref = [q0; w0];
options = odeset('AbsTol', tol);
[timeSol, zC] = ode45(odeFunc, [0 T], z0Ref, options);
[err, times, sol] = variableStepComparison(vecField, action, z0, T, tol) ;
zC = zC';
zC = reorder(zC);
plot(TT(1 : end - 1), diff(TT), 'k--', 'linewidth', 2);
hold on;
plot(timeSol(1 : end - 1), diff(timeSol), 'b:', 'linewidth', 2);
hold on;
plot(timeSol(1 : end - 1), 0 * timeSol(1 : end - 1) + dt, 'r-', 'linewidth', 2);

legend('RKMK(5, 4)', 'ODE45', 'RKMK5', 'FontSize', 30, 'Location', 'northwest');
xlim([0 T])
xlabel("Time", 'fontsize', 30);
ylabel("Stepsize", 'fontsize', 30);
ax = gca;
ax.XAxis.FontSize = 30;
ax.YAxis.FontSize = 30;
stringa = 'Comparison of stepsize variation with N = '+string(P)+' pendulums';
