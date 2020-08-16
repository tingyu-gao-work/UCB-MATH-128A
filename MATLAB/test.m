interval = linspace(0, 2 * pi, 100); 
[X,Y] = meshgrid(interval, interval);

% f = @(x,y) (2 * sin(x) * cos(y));
% Z = arrayfun(f,X,Y);
Z = 2 * sin(X) .* cos(Y);

f = figure;
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
zlabel('$z$', 'Interpreter', 'latex');

hold on
surf(X,Y,Z);
view(3)
hold off

lg = legend('$2 \sin{(x)} \cos{(y)}$');
lg.Interpreter = 'latex';