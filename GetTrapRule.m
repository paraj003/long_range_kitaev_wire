function [x, w] = GetTrapRule(N, t0, tf)
    % Computes the points and weights for the Trapazoidal quadrature
    % rule from t0 to tf with N function evaluations.
    % Should converge >~ N^-2.
    N = N - 1; delta = (tf - t0) / N;
    x = zeros(1, N+1); w = zeros(1, N+1);
    x(1) = t0; x(end) = tf;
    w(1) = delta / 2.0; w(end) = delta / 2.0;
    for n=1:N-1
        x(n+1) = t0 + n * delta;
        w(n+1) = delta;
    end
end
