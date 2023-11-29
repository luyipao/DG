x=1
f = @(x) x;
g = @(t) t;
h = @(x) f(x) + g(x);
h(1)