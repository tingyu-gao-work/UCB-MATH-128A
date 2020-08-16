t = 0:1:5;
x = [1, 1.5, 2, 2, 2.5, 2.5];
y = [1, 0.5, 1, 1.5, 1.5, 1];

[bx, cx, dx] = ncspline(t, x);
[by, cy, dy] = ncspline(t, y);

fprintf('\nX:');
fprintf('\na:'); disp(x(1:length(x)-1))
fprintf('\nb:'); disp(bx)
fprintf('\nc:'); disp(cx)
fprintf('\nd:'); disp(dx)

fprintf('\nY:');
fprintf('\na:'); disp(y(1:length(y)-1))
fprintf('\nb:'); disp(by)
fprintf('\nc:'); disp(cy)
fprintf('\nd:'); disp(dy)

T = linspace(0,5,1000);
xx = splineeval(t, x(1:length(x)-1), bx, cx, dx, T);
yy = splineeval(t, y(1:length(y)-1), by, cy, dy, T);

plot(xx, yy, x, y, 'o');
axis equal
grid on