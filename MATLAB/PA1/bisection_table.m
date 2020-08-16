function p = bisection_table(f, a, b, tol)
% Solve f(p) = 0 using the bisection method.
% Like bisection, but printing a convergence table.

% Print header
fprintf(' n         a             b             p           f(p)    \n');
fprintf('-----------------------------------------------------------\n');

for n = 1:100 % Terminate after 100 iterations (even if not converged)
    p = (a+b) / 2;
    if p-a < tol, break; end
    fprintf('%2d  %12.8f  %12.8f  %12.8f  %12.8f\n', n, a, b, p, f(p));
    if f(a)*f(p) > 0
        a = p;
    else
        b = p;
    end
end
