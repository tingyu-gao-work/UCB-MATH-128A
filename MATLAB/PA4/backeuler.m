function [t,w] = backeuler(f, dfdy, a, b, alpha, N, maxiter, tol)
h = (b-a)/N;
t = linspace(a, b, N+1);

w = zeros(1, N+1);
w(1) = alpha;

for i = 1:N
    w(i+1) = newton(f, dfdy, h, t(i+1), w(i), maxiter, tol, i);
end
