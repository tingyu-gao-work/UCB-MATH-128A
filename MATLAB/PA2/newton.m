function p = newton(f, df, p0, tol)
% Solve f(p) = 0 using Newton's method.

while 1
    p = p0 - f(p0)/df(p0);
    if abs(p-p0) < tol, break; end
    p0 = p;
end