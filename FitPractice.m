x = [1,4,9];
y = sqrt(x);
p = polyfit(x,y,2);
%disp(polyval(p,3))

inter = linspace(0,10,100);

f = figure;
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');

hold on
p1 = plot(inter, sqrt(inter));
p2 = plot(inter, polyval(p,inter));
hold off

grid on
lg = legend('$y = \sqrt{x}$','$y = p(x)$');
lg.Interpreter = 'latex';