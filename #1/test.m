% [PN, DPN] = Legendre(0, -1, 1);
phi_jj = @(x) 1;
phi_m = @(x) 2;
quadgk(@(x)  arrayfun(@(x) phi_jj(x) .* phi_m(x),x),1, 2)

function [phi, diff_phi] = base_func(k,xa,xb)
h = xb - xa;
phi = @(x) ((x-xa).^k / (h^k)) * (k>=1) + 1.0 * (k ==0);
% derivation of base function phi
diff_phi = @(x) (k >=1) * (k * (x-xa).^(k-1) / h^k);
end

%生成区间[a,b]上的标准Legendre多项式
function [Pn,DPn] = L(n,a,b)
Q = @(x) 0.*x;
DQ = @(x) 0.*x;
R = @(x) 0.*x + 1;
DR = @(x) 0.*x;
P = @(x) ((2*0 + 1).*x.*R(x) - 0.*Q(x))/(0 + 1);
DP = @(x) ((2*0 + 1).*(x.*DR(x) + R(x)) - 0.*DQ(x))/(0 + 1);
if n == -1
    Pn = Q;
    DPn = DQ;
elseif n == 0
    Pn = R;
    DPn = DR;
elseif n == 1
    Pn = P;
    DPn = DP;
else
    for i = 3:n + 1
        k = i - 2;
        Q = R;
        DQ = DR;
        R = P;
        DR = DP;
        P = @(x) ((2*k + 1).*x.*R(x) - k.*Q(x))/(k + 1);
        DP = @(x) ((2*k + 1).*(x.*DR(x) + R(x)) - k.*DQ(x))/(k + 1);
    end
    Pn = P;
    DPn = DP;
end
%尺度变换
c = (a + b)/2;
h = (b - a)/2;
Pn = @(x) sqrt((2*n + 1)/(2*h))*Pn((x - c)/h);
DPn = @(x) sqrt((2*n + 1)/(2*h^3))*DPn((x - c)/h);
end