function p = newton_plot(f, df, p0, xmin, xmax)
% Solve f(p) = 0 using Newton's method.
% Like newton, but no convergence check (always 5 iterations)
% and plotting the graph and the iterations.

% Print header
fprintf(' n         p          |p-p0|   \n');
fprintf('-------------------------------\n');

% Plot f(x)
xx = linspace(xmin, xmax, 1000);
plot(xx, f(xx), 'linewidth',2)
line([xmin,xmax], [0,0], 'linewidth',2, 'color','b');
grid on

for n = 1:5
    p = p0 - f(p0)/df(p0);
    fprintf('%2d  %12.8f  %12.8f\n', n, p, abs(p-p0));
    line([p0,p0], [0,f(p0)], 'color','r', 'linewidth',2);
    line(p0, f(p0), 'marker','.', 'color','r'); pause(1);
    line([p0,p], [f(p0),0], 'color','r', 'linewidth',2); pause(1);
    p0 = p;
end