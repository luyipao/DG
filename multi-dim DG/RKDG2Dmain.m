
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


xa = 0;
xb = 2*pi;
ya = 0;
yb = 2*pi;
tb = 3;
k = 2;
N = 40;
CFL = 0.05;
[Q,C] = RKDG2D(xa,xb,ya,yb,tb,k,N,CFL);

function [Q, C] = RKDG2D(xa, xb, ya, yb, tb, k, N, CFL) 
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

classicalLegendreBasis2DLeft = zeros(6,1);
for n = 1:6
    [i,j] = findCoordinate(n);
    classicalLegendreBasis2DLeft(n) = classicalLegendreBasis2D{n}(-1,-1);

end

hx = (xb - xa) / N;
hy = (yb - ya) / N;
X = linspace(xa, xb, N + 1);
Y = linspace(ya, yb, N + 1);
X = X(:);
Y = Y(:);
T = 0;
dim = (k + 1) * (k + 2) / 2;
A = zeros(dim);
B = zeros(dim);
D = zeros(dim);
E = zeros(dim);

classicalLegendreBasis2DLeft = zeros(6,1);
for n = 1:6
    [i,j] = findCoordinate(n);
    classicalLegendreBasis2DLeft(n) = sqrt(2*i+1) * sqrt(2*j+1) / sqrt(hx*hy) * classicalLegendreBasis2D{n}(-1,-1);
end

for m = 1:dim
    [mi, mj] = findCoordinate(m);
    mi = mi - 1;
    mj = mj - 1;
    for n = 1:dim
        [ni, nj] = findCoordinate(n);
        ni = ni - 1;
        nj = nj - 1;
        A(m, n) = 4 / ((2 * mi + 1) * (2 * mj + 1)) * (mi == ni) * (mj == nj);
        B(m, n) = -sqrt((2 * mj + 1) * (2 * nj + 1)) / hy * (mi == ni);
        B(m, n) = B(m, n) + sqrt((2 * mi + 1) * (2 * ni + 1)) / hx * (mj == nj);
        B(m, n) = B(m, n) - (ni > mi) * (mod(ni - mi, 2) == 1) * 2 * (mj == nj) ...
            - (nj > mj) * (mod(nj - mj, 2) == 1) * 2 * (mi == ni);
        D(m, n) = (-1)^(mi + ni + 1) * (mj == nj);
        E(m, n) = (-1)^(mj + nj) * (mi == ni);
    end
end

% 初始化 C 为数值数组
C = zeros(dim, (N + 1) ,(N + 1)); % 这里将 C 设为一个二维数组
Q = zeros(N, N, 10000);
for i = 2:N+1
    for j = 2:N+1
        c = zeros(dim, 1);
        for kk = 1:dim
            c(kk) = quad2d(@(x, y) classicalLegendreBasis2D{kk}(x, y) .* f(x, y), X(i), X(i + 1), Y(j), Y(j + 1));
        end
        C(:, i , j) = c; % 将 c 存储到 C 的相应位置
        Q(i,j,1) = dot(c,classicalLegendreBasis2DLeft(1:dim));
    end
end


Qtrue = zeros(N,N);
for i = 1:N
    for j = 1:N
Qtrue(i,j) = f(X(i),Y(j));
    end
end

figure(1)
surf(X(1:N),Y(1:N),Q(1:40,1:40,1))
figure(2)
surf(X(1:N),Y(1:N),Qtrue)

% 周期边界条件
C(:,N+1,:) = C(:,1,:);
C(:,:,N+1) = C(:,:,1);

t = CFL * (sqrt(hx * hy))^((k + 1) / 3);
while (T < tb) 
    C = GK3(C, A, B, D, E, t);
    Q = zeros(N, N);
    for i = 2:N + 1
        for j = 2:N + 1
            for kk = i:dim
                [m, n] = findCoordinate(kk);
                Q(i - 1, j - 1) = Q(i - 1, j - 1) + C(kk, i, j) * classicalLegendreBasis2D{kk}(-1, -1) * (-1)^(m + n) * sqrt(2 * m + 1) * sqrt(2 * n + 1) / sqrt(hx * hy);
            end
        end
    end
end
end

function C = GK3(C, A, B, D, E, t)
k1 = C + t * L(C, A, B, D, E);
k2 = (3 / 4) * C + (1 / 4) * k1 + (1 / 4) * t * L(k1, A, B, D, E);
C = (1 / 3) * C + (2 / 3) * k2 + (2 / 3) * t * L(k2, A, B, D, E);
end

function C = L(C, A, B, D, E)
tempA = inv(A);
[~,~,N] = size(C);
N = N - 1;
for i = 2:N + 1
    for j = 2:N + 1
        C(:, i, j) = -tempA * B * C(:, i, j) - tempA * D * C(:, (i - 1) * (N + 1) + j - 1) - tempA * E * C(:, (i - 2) * (N + 1) + j);
    end
end
C(:,1,:) = C(:,end,:);
C(:,:,1) = C(:,:,end);
end

function y = f(x, y)
    y = sin(x) .* sin(y);
end

function [i, j] = findCoordinate(n)
    k = (sqrt(8 * n + 1)) - 1;
    k = ceil(k / 2);
    n = n - k * (k + 1) / 2;
    i = 1 - n; 
    j = k + n; 
end
