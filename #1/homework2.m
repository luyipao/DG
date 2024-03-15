f = @(x) 2 * pi * cos(2* pi *x);
% u(x,t) = sin(2pi(x - t))
%input: 区间[xa,xb]；cell宽度size；解空间多项式最高阶degree；t=0时函数u(x,0)=u_0
% output: 系数矩阵C，其列向量c_i表示在cell I_i上的u的系数
u_0 =@(t)  sin(2*pi * t);
u = @(x,t) sin(2*pi *(x-t));
degree = 0;
kt = 1;
xa = 0;
xb = 1;
ta = 0;
tb = 1;
CFL = 0.1;
n = 3;
h = geospace(1/10,1/2,n);
for ii = 1:n
    [u_hh,U_h] = RKDG(f,xa,xb,ta,tb,u_0,h(ii),degree,kt,CFL);
    E(ii) = analysisDG(u,xa,xb,ta,tb,h(ii),degree,kt, U_h,CFL);
end
KL1 = E(1:end-1) ./ E(2:end);
KL1 = log(KL1) / log(2);

function L1 = analysisDG(u,xa,xb,ta,tb,h,k,kt,U_h,CFL)
X = xa:h:xb;
h_t = CFL * h ;
T = ta:h_t:tb;
U = zeros(length(T), length(X));
for ii = 1:length(X)
    for tt = 1:length(T)
        U(tt,ii) = u(X(ii),T(tt));
    end
end
E_h = U - U_h;
L1 = norm(E_h,1) / (length(X) * length(T));
disp(['h is ', num2str(h), '; k is: ', num2str(k), '; kt is: ', num2str(kt), '; error is ', num2str(sum(sum(abs(E_h))) / (length(X)*length(T)))]);
end
%input: h 网格宽度；kt RK法的阶；degree 多项式最高阶；CFL cfl数
% output 时间t_1到t_n的函数u_hh；网格上的数值解矩阵D
function [u_hh,D] = RKDG(f,xa,xb,ta,tb,u_0,h,degree,kt,CFL)
    X = xa:h:xb;
    h_t = CFL * h;
    T = ta:h_t:tb;
    N_t = length(T) - 1;
    u_hh = cell(length(T),1);% t_n 上近似[xa,xb]的u_h
    [~,C] = DG(f,xa,xb,u_0(xa),h,degree);
    for tt = 1:N_t
        u_hh{tt} = builder(C,X);
        if kt == 1
            [~,C] = Euler_solve(u_hh{tt}, h_t, X, C);
        elseif kt == 2
            [~,C] = RK2(u_hh{tt}, h_t, X, C);
        elseif kt == 3
            [~,C] = RK3(u_hh{tt}, h_t, X, C);
        end
    end
    u_hh{length(T)} = builder(C,X);
    D = zeros(length(T),length(X));
    for tt = 1:length(T)
        D(tt,:) = u_hh{tt}(X);
    end
end


%input: 区间[xa,xb], 基函数维数+1  K；算子L；时间间隔t；函数u；步长 h
%output:  欧拉迭代后的系数矩阵C；阶段函数 u_hh；
function [u_hh, C] = Euler_solve(u, delta_t, X, C)
[K, N] = size(C);
A = zeros(K,K);
for ii = 1:N
    for m = 1:K
        for jj = 1:K
            [phi_m,~] = base_func(m-1,X(ii),X(ii+1));
            [phi_jj,~] = base_func(jj-1,X(ii),X(ii+1));
            A(m,jj) = quadgk(@(x) phi_m(x) .* phi_jj(x), X(ii), X(ii+1));
        end
    end
    b = L(u, X(ii), X(ii+1), K);
    C(:,ii) = C(:,ii) + delta_t  * (A\b);
end
u_hh = builder(C,X);
end

% 
function [u_hh, C] = RK2(u_h,delta_t,X,C)
[u1,C1] = Euler_solve(u_h, delta_t,X,C);
[u_temp, C_temp] = Euler_solve(u1, delta_t,X,C1);
C = 0.5 * (C + C_temp);
u_hh = @(t) 0.5 * u_h(t) + 0.5 * u_temp(t);
end

function [u_hh, C] = RK3(u_h,  delta_t,X,C)
[u1,C1] = Euler_solve(u_h,delta_t,X, C);
[~,C_temp] = Euler_solve(u1,delta_t,X,C1);
C2 = 3 / 4 * C + 1 /4 * C_temp;
u2 = builder(C2,X);
[u_temp,C_temp] = Euler_solve(u2,delta_t, X, C2);
C = 1 / 3 * C + 2/3 * C_temp;
u_hh = @(t) 1/3* u_h(t) + 2/3 * u_temp(t);
end

%函数u_hh
function u_hh = builder(C,X)
[K, N] = size(C);
u_h = repmat({@(t) 0}, 1, N);
u_hh = @(t) 0;
for ii = 1:N
    for jj = 1:K
        [phi,~] = base_func(jj-1,X(ii),X(ii+1));
        u_h{ii} = @(t) u_h{ii}(t) + C(jj,ii) * phi(t);
    end
    u_hh  = @(t) u_hh(t) + ((X(ii) < t) & (t <= X(ii+1))) .* u_h{ii}(t);
end
u_hh = @(t) u_hh(t) + (t==X(1)) .* u_h{1}(t);
end


% input: 多项式的阶k, cell [a,b]
% output: 基函数phi和其导数diff_phi
function [phi, diff_phi] = base_func(k,xa,xb)
h = xb - xa;
phi = @(x) ((x-xa).^k / (h^k)) * (k>=1) + 1.0 * (k ==0);
% derivation of base function phi
diff_phi = @(x) (k >=1) * (k * (x-xa).^(k-1) / h^k);
end

% u函数；xa cell左端点；xb cell右端点；K 基函数的维数
%output: b RK法中的f(t,y)
function b = L(u,xa,xb,K)
b = zeros(K,1);
for ii = 1:K
    [phi,diff_phi] = base_func(ii-1,xa,xb);
    b(ii) = quadgk(@(x) u(x) .* diff_phi(x), xa, xb) - u(xb) * phi(xb) + u(xa) * phi(xa);
end
end

function h = geospace(a,r,n)
h = zeros(1, n);  % 初始化等比数列
h(1) = a;  % 设置起始值a
for i = 2:n
    h(i) = h(i-1) *r;  % 计算下一个元素，公比为 r
end
end