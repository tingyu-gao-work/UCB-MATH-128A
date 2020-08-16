t = 0:1:5;
x = [1, 1.5, 2, 2, 2.5, 2.5];
y = [1, 0.5, 1, 1.5, 1.5, 1];

[bx, cx, dx] = ncspline(t, x);
[by, cy, dy] = ncspline(t, y);

y_t_0 = @(s) (splineeval(t, y(1:length(y)-1), by, cy, dy, s) - 1.2);
dy_t_0 = @(s) diffsplineeval(t, y(1:length(y)-1), by, cy, dy, s);
dx_t_0 = @(s) diffsplineeval(t, x(1:length(x)-1), bx, cx, dx, s);
int_func = @(s) sqrt(dx_t_0(s).^2 + dy_t_0(s).^2);

t1 = newton(y_t_0, dy_t_0, 2, 1e-8);
t2 = newton(y_t_0, dy_t_0, 5, 1e-8);

n = [16,32,64,128,10000];
L = zeros(1,5);
for i = 1:5
    L(i) = Composite_Trapezoidal(int_func,t1,t2,n(i));
end
fprintf('[L_16, L_32, L_64, L_128, L]\n')
fprintf('%12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n', L(1), L(2), L(3), L(4), L(5));

Delta_L_abs = abs(L(1:length(L)-1) - L(length(L)));
h = (t2-t1)./n(1:length(L)-1);
g = fittype('a*x^b');
f = fit(h(:), Delta_L_abs(:), g, 'StartPoint', [1, 2]);
disp(f);


hold on
plot(f);
scatter(h, Delta_L_abs, 'filled');
hold off

set(gca,'xscale','log', 'yscale','log');
xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');

grid on

lg = legend('Fitted Curve','$|L_n - L|$','Location','southeast');
lg.Interpreter = 'latex';
function p = Composite_Trapezoidal(f,a,b,n)
h = (b-a)/n;
j = 1:1:n-1;
x = a * ones(1,n-1) + h * j;
p = h/2 * (f(a) + f(b)) + h * sum(f(x));
end