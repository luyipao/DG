clf
%% initial Function
x = linspace(0, 1, 1000); % test interval [0, 1]
xp = rand(10,10); % test matrix
% draew
hold on
plot(x,initialFunction(x));
scatter(xp,initialFunction(xp));
legend("initial function");
hold off

