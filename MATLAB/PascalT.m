T = {};
for i = 1:20
    T{end + 1} = calculate(i);
end
disp(T{10})
disp(T{20})

function p = calculate(n)
    if n == 1
        p = [1];
    elseif n == 2
        p = [1,1];
    else
        row = calculate(n-1);
        p = [row, 0] + [0, row];
    end
end