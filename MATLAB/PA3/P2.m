f = @(t, y) fpend(y);

[t1, w1] = rk4(f, 0, 100, [1,1,0,0], 2000);
[t2, w2] = rk4(f, 0, 100, [pi,0,0,1e-10], 2000);
[t3, w3] = rk4(f, 0, 100, [2,2,0,0], 2000);
[t4, w4] = rk4(f, 0, 100, [2,2 + 1e-3,0,0], 2000);


hold on
plot(t1, w1(:,2));
plot(t2, w2(:,2));
plot(t3, w3(:,2));
plot(t4, w4(:,2));
hold off

xlabel('$t$', 'Interpreter','latex');
ylabel('$\theta_2$', 'Interpreter','latex');

grid on
lg = legend('Case1','Case2','Case3','Case4');