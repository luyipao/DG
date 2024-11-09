% 定义经典勒让德基函数
classicalLegendreBasis = {
    @(x) 1 + 0 * x, ... % P_0(x)
    @(x) x, ...         % P_1(x)
    @(x) 1/2 * (3 * x.^2 - 1), ... % P_2(x)
    @(x) 1/2 * (5 * x.^3 - 3 * x), ... % P_3(x)
    @(x) 1/8 * (35 * x.^4 - 30 * x.^2 + 3) % P_4(x)
};

% 定义对应的导函数
classicalLegendreDerivatives = {
    @(x) 0 * x, ... % P_0'(x)
    @(x) 1 + 0 * x, ... % P_1'(x)
    @(x) 3 * x, ... % P_2'(x)
    @(x) 3/2 * (5 * x.^2 - 1), ... % P_3'(x)
    @(x) 1/2 * (140 * x.^3 - 60 * x) % P_4'(x)
};



classicalLegendreBasis2D = {
    @(x, y) classicalLegendreBasis{1}(x) .* classicalLegendreBasis{1}(y), ... % 1
    @(x, y) classicalLegendreBasis{2}(x) .* classicalLegendreBasis{1}(y), ... % x
    @(x, y) classicalLegendreBasis{1}(x) .* classicalLegendreBasis{2}(y), ... % y
    @(x, y) classicalLegendreBasis{3}(x) .* classicalLegendreBasis{1}(y), ... % x^2
    @(x, y) classicalLegendreBasis{2}(x) .* classicalLegendreBasis{2}(y), ... % xy
    @(x, y) classicalLegendreBasis{1}(x) .* classicalLegendreBasis{3}(y)      % y^2
};

function [i, j] = findCoordinate(n)
    % basis: 基函数排序为P_0(x)P_0(y), P_1(x)P_0(y),
    % P_0(x)P_1(y),P_2(x)P_0(y),P_1(x)P_1(y),P_0(x)P_2(y)
    % input: n 从1开始，表示basis的第n个基函数
    % output: [i,j] 从[1,1]开始表示基函数P_{i-1}(x)P_{j-1}(y)
    % k 是n所对应的最大次数也就是 i + j = k + 1.
    k = (sqrt(8*n+1)) - 1;
    k = ceil(k / 2);
    n = n - k * (k + 1) / 2;
    i = 1 - n; 
    j = k + n; 
end



% for n = 1:6
%     [i,j] = findCoordinate(n);
%     disp(calssicalLegendreBasis{n});
% end
% 示例：计算在 x = 0.5 处的导函数值
% x = 0.5;
% derivativeValues = cellfun(@(f) f(x), classicalLegendreDerivatives);
% disp(derivativeValues);

% 示例：计算在 x = 0.5 处的基函数值
% x = 0.5;
% values = cellfun(@(f) f(x), classicalLegendreBasis);
% disp(values);
