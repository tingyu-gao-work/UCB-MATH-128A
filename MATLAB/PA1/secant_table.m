function p = secant_table(f, p1, p2, tol)
% Solve f(p) = 0 using secant method.

% Print header
fprintf(' n        p       |p_n-p_{n-1}|        f(p)    \n');
fprintf('-----------------------------------------------\n');

n = 0;
while 1
    n = n + 1;
    p3 = p2 - f(p2)/((f(p2) - f(p1))/(p2-p1));
    delta = abs(p3-p2);
    fprintf('%2d  %12.8f  %12.8e  %12.8e\n', n, p3, delta, f(p3));
    if delta < tol || abs(f(p3))<tol
        p = p3;
        break
    end
    p1 = p2;
    p2 = p3;
end