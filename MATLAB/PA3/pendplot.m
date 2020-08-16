function pendplot(th1,th2)
% Plot a double pendulum with angles th1,th2

% x,y position of first weight
x1 = sin(th1);
y1 = -cos(th1);

% x,y position of second weight
x2 = x1 + sin(th2);
y2 = y1 - cos(th2);

% Plot using lines and circles
cla, axis equal, axis([-2.1,2.1,-2.1,2.1]), grid on
line([-.2,.2], [0,0], 'linewidth',5, 'color','b')
line([0,x1,x2], [0,y1,y2], 'linewidth',3, 'color','k')
line([x1,x2], [y1,y2], 'linestyle','none', 'marker','.', ...
     'markersize',32, 'color','r')
pause(0.01)