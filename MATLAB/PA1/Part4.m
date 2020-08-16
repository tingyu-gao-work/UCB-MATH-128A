syms f_1(x) f_2(x)

f_1(x) = sin(x) - exp(-x);
f_2(x) = (sin(x^2)) / (10 + x^2) - 1/50 * exp(-x/10);
df_1 = diff(f_1,x);
df_2 = diff(f_2,x);

f_1 = matlabFunction(f_1);
f_2 = matlabFunction(f_2);
df_1 = matlabFunction(df_1);
df_2 = matlabFunction(df_2);




f_1_zeros = findmanyzeros(f_1, 0, 10, 50, 1e-10);
f_2_zeros = findmanyzeros(f_2, 0, 10, 50, 1e-10);
df_1_zeros = findmanyzeros(df_1, 0, 10, 50, 1e-10);
df_2_zeros = findmanyzeros(df_2, 0, 10, 50, 1e-10);


t = linspace(0,10);
fig1 = figure;
xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');

hold on
plot(t, f_1(t));
plot(t, zeros(1,length(t)))
scatter(f_1_zeros, zeros(1, length(f_1_zeros)))
scatter(df_1_zeros, f_1(df_1_zeros), '^')
hold off

grid on
lg = legend('$f_1(x)$','$y = 0$','Zeros','Extrema','Location','southeast');
lg.Interpreter = 'latex';




t = linspace(0,10,500);
fig2 = figure;
xlabel('$x$', 'Interpreter','latex');
ylabel('$y$', 'Interpreter','latex');

hold on
plot(t, f_2(t));
plot(t, zeros(1,length(t)))
scatter(f_2_zeros, zeros(1, length(f_2_zeros)))
scatter(df_2_zeros, f_2(df_2_zeros), '^')
hold off

grid on
lg = legend('$f_2(x)$','$y = 0$','Zeros','Extrema');
lg.Interpreter = 'latex';