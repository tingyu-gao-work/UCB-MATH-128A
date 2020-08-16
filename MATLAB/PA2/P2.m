t = 0:1:5;
y = [1, 0.5, 1, 1.5, 1.5, 1];
[by, cy, dy] = ncspline(t, y);

y_t_0 = @(s) (splineeval(t, y(1:length(y)-1), by, cy, dy, s) - 1.2);
dy_t_0 = @(s) diffsplineeval(t, y(1:length(y)-1), by, cy, dy, s);

p1 = newton(y_t_0, dy_t_0, 2, 1e-8);
p2 = newton(y_t_0, dy_t_0, 5, 1e-8);

fprintf('t1 = %.8f \n', p1);
fprintf('t2 = %.8f \n', p2);