function output = checkCondition(i, j)
    output = (j > i) * (mod(j - i, 2) == 1) * 2; % 当 j = i + 1, i + 3, i + 5, ... 时输出 2，否则输出 0
end
for i = 1:6
    for j = 1:6
        A(i,j) = checkCondition(i, j);
    end
end
