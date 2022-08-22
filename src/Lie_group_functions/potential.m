function val = potential(q,L,m)

    P = length(L);
    g = 9.81;
    e3 = [0;0;1];
    
    val = 0;
    for i = 1:P
        val = val + L(i)*g*sum(m(i:end))*e3'*q(3*i-2:3*i);
    end

end