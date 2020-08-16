function p = newton_table(f, df, p0, niter)
% Solve f(p) = 0 using Newton's method.
% Like newton, but no convergence check (always niter iterations).

% Print header
fprintf(' n         p          |p-p0|   \n');
fprintf('-------------------------------\n');

for n = 1:niter
    p = p0 - f(p0)/df(p0);
    fprintf('%2d  %12.8f  %12.8f\n', n, p, abs(p-p0));
    p0 = p;
end