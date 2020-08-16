theta = linspace(0, 2 * pi, 1000);
x = 16 * (sin(theta)).^3;
y = 13 * cos(theta) - 5 * cos(2 * theta) - 2 * cos(3 * theta) - cos(4 *theta);


f = figure;
xlabel('$\theta$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');

hold on
p1 = plot(theta, x);
p2 = plot(theta, y);
hold off

grid on

lg = legend('$f(x)$', '$g(x)$');
lg.Interpreter = 'latex';
