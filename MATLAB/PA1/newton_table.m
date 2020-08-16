function p = newton_table(f, df, p0, tol)
% Solve f(p) = 0 using Newton's method.

% Print header
fprintf(' n        p           |p-p0|           f(p)    \n');
fprintf('-----------------------------------------------\n');

n = 0;
while 1
    n = n + 1;
    p = p0 - f(p0)/df(p0);
    delta = abs(p-p0);
    fprintf('%2d  %12.8f  %12.8e  %12.8e\n', n, p, delta, f(p));
    if delta < tol %|| abs(f(p))<tol
        break
    end
    p0 = p;
end