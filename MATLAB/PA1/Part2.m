f = @(x) cos(x) - x;
%df = @(x) -sin(x) -1;
findzero(f, 0, 2, 1e-15)
%newton_table(f,df,0,1e-10)
secant_table(f, 2, 0, 1e-15)