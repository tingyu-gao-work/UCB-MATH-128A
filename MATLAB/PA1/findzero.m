function p = findzero(f, a, b, tol)
% Solve f(p) = 0 using the bisection and secant method.


% Printing a convergence table.
% Print header
%fprintf(' n         a             b             p         abs(b-a)  \n');
%fprintf(' n         a             b             p           f(p)    \n');
fprintf(' n         p         abs(b-a)           f(p)    \n');
fprintf('------------------------------------------------\n');
%fprintf('-----------------------------------------------------------\n');


w = 1;
for n = 1:100
    p = a + (w * f(a) * (a - b)) / (f(b) - w * f(a));
    %fprintf('%2d  %12.8f  %12.8f  %12.8f  %12.8e\n', n, a, b, p, abs(b-a));
    %fprintf('%2d  %12.8f  %12.8f  %12.8f  %12.8e\n', n, a, b, p, f(p));
    fprintf('%2d  %12.8f  %12.8e  %12.8e\n', n, p, abs(b-a), f(p));
    
    if f(p) * f(b) > 0
        w = 1/2;
    else
        w = 1;
        a = b;
    end
    
    b = p;
    if abs(b-a) < tol || abs(f(p)) < tol
        break
    end
end

end