clf
%% initial Function
x = linspace(0, 1, 1000); % test interval [0, 1]
xp = rand(10,10); % test matrix
% draew
hold on
plot(x,N_D(x));
scatter(xp,N_D(xp));
legend("initial function");
hold off
%% basisPolys
clf
%  系数是正交归一化的系数，基目前必须是勒让德基函数
basisFuncs = {@(x) 1 + 0 * x, @(x) x, @(x) 0.5 * (3 * x.^2 + 1)};
mesh = [-1, 1];
interval = [-1, 1];
coeffs = [1.4142, 0.8165, 0.6325]';
degree = 2;
f = basisPolys(mesh,coeffs,degree,basisFuncs,interval);
X = linspace(-1, 1);
hold on
plot(X, 1 + X + 0.5 * (3 * X.^2 + 1));
plot(X, f.solve(X), '--');
hold off
