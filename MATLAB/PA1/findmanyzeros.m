function p = findmanyzeros(f, a, b, n, tol)
Zeros = [];

x = linspace(a, b , n+1);
for i = 1:n
    if f(x(i)) * f(x(i+1)) < 0
        temp = findzero(f, x(i), x(i+1), tol);
        Zeros(length(Zeros) + 1) = temp;
    end
end
p = Zeros;
end