function m = create_m(N, B, J, eta)
    % create the matrix m where H_eta = A^dag m A (+ const) for
    % nearest neighbor interaction.
    % eta = 0: open boundary conditions.
    % eta = 1: periodic boundary conditions.
    % eta = -1: antiperiodic boundary conditions.
    % see appendix C.1.

    m = zeros(2*N, 2*N); a = J / 2;
    for i=1:N
        m(2*i-1, 2*i-1) = B / 2.0;
        m(2*i, 2*i) = -B / 2.0;
        j = i + 1;
        if j <= N
            [m(2*i-1, 2*j-1), m(2*j-1, 2*i-1)] = deal(-a);
            [m(2*i-1, 2*j), m(2*j, 2*i-1)] = deal(-a);
            [m(2*i, 2*j-1), m(2*j-1, 2*i)] = deal(a);
            [m(2*i, 2*j), m(2*j, 2*i)] = deal(a);
        end
    end
    i = N; j = 1;
    [m(2*i-1, 2*j-1), m(2*j-1, 2*i-1)] = deal(-eta*a);
    [m(2*i-1, 2*j), m(2*j, 2*i-1)] = deal(-eta*a);
    [m(2*i, 2*j-1), m(2*j-1, 2*i)] = deal(eta*a);
    [m(2*i, 2*j), m(2*j, 2*i)] = deal(eta*a);
end

