function p = newton(f, dfdy, h, t, w, maxiter, tol, i)
% Solve w_{i+1} = 0 using Newton's method.
% t = t_{i+1}; w = w_i;

g = @(x) h * f(t, x) + w - x;
dg = @(x) h * dfdy(t, x) - 1;
p0 = w;
iter = 0;

% Print header
fprintf('This is the table for step %d.\n', i);
fprintf(' n         p          |p-p0|   \n');
fprintf('-------------------------------\n');

while 1
    iter = iter + 1;
    p = p0 - g(p0)/dg(p0);
    delta = abs(p-p0);
    
    fprintf('%2d  %12.8f  %12.8f\n', iter, p, delta);
    if delta < tol || iter == maxiter
        fprintf('\n');
        break; 
    end
    p0 = p;
end