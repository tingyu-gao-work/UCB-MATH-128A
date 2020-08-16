N = [2000, 4000, 8000, 16000];
h = 100./N;
error = zeros(1,4);

f = @(t, y) fpend(y);
[t0, w0] = rk4(f, 0, 100, [1,1,0,0], 100000);
t2_0 = w0(:,2);
real = t2_0(end);
fprintf('h = %.5f: The exact value is roughly %.8f.\n', 0.001, real);

for i = 1:4
    [t, w] = rk4(f, 0, 100, [1,1,0,0], N(i));
    t2 = w(:,2);
    est = t2(end);
    fprintf('h = %.5f: The value is roughly %.8f.\n', 100/N(i), est);
    error(i) = abs(est - real);
end

P = polyfit(log(h), log(error), 1);
f = @(x) exp(P(2))*x.^P(1);
xx = linspace(0.001, 0.1, 100);
fprintf('The slope is roughly %.2f.\n', P(1))

hold on
plot(xx, f(xx));
scatter(h, error, 'filled');
hold off

set(gca,'xscale','log', 'yscale','log');
xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');
grid on

lg = legend('Fitted Curve','$|w_f - \theta_2(100)|$', 'Location','southeast');
lg.Interpreter = 'latex';