dudt = @(t, u) f(t,u);
[t, w] = rk4(dudt, 0, 1/2, [1/3,1/3], 20);

tt = linspace(0, 1/2);
u1 = @(t) 2/3 * t + 2/3 * exp(-t) - 1/3 * exp(-100 * t);
u2 = @(t) -1/3 * t - 1/3 * exp(-t) + 2/3 * exp(-100 * t);
y1 = u1(tt);
y2 = u2(tt);

hold on
plot(tt, y1)
plot(t, w(:,1))
hold off


function du = f(t,u)
    u1 = u(1);
    u2 = u(2);
    
    du1_dt = 32 * u1 + 66 * u2 + 2/3 * t + 2/3;
    du2_dt = -66 * u1 - 133 * u2 - 1/3 * t - 1/3;
    
    du = [du1_dt, du2_dt];
end