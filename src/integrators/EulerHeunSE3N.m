function sol = EulerHeunSE3N(vecField,action,p,h)

    %Computes one time step update with Lie Euler Heun's method

    k0 = zeros(length(p),1);
    k1 = vecField(k0,p);
    k2 = vecField(h*k1, p);
    sol = action(exponentialSE3N(h/2*(k1+k2)),p);

end