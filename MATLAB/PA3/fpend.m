function ydot = fpend(y)
    t1 = y(1);
    t2 = y(2);
    w1 = y(3);
    w2 = y(4);

    dt1_dt = w1;
    dt2_dt = w2;
    dw1_dt = (-3*sin(t1) - sin(t1 - 2*t2)...
        - 2*sin(t1-t2) .* (w2.^2 + w1.^2 .* cos(t1-t2)))...
        ./...
        (3 - cos(2*(t1-t2)));
    dw2_dt = (2*sin(t1-t2) .* (2*w1.^2 + 2*cos(t1) + w2.^2 .* cos(t1-t2)))...
        ./...
        (3 - cos(2*(t1-t2)));
    ydot = [dt1_dt, dt2_dt, dw1_dt, dw2_dt];
end
