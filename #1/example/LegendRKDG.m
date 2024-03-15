%间断有限元方法求解一维双曲守恒率：
%u_t + u_x = 0,初值条件为 u(x,0) = f(x),
%边界条件取周期边界，时间上使用龙格库塔方法
function [X,T,Q,C] = LegendRKDG(xa,xb,tb,k,kt,N,CFL)
%输入：空间求解域[xa,xb]，时间求解域[0,tb]，
%空间剖分段数N(时间步长自适应生成)，
%分片多项式的最高次数k，基函数使用尺度变换后的标准Legendre多项式
%时间阶数kt，可选用3,4阶龙格库塔方法，以及CFL数
%输出：空间网格X，时间网格T，近似解矩阵Q，
%以及基函数的系数矩阵C

%初始化
h = (xb - xa)/N;
X = xa:h:xb;
X = X';
T = 0;

%构造刚度矩阵(由于选取了标准正交基，质量矩阵是单位阵)

G = zeros(k + 1,k + 1);
H = zeros(k + 1,k + 1);
for i = 1:k + 1
    for j = 1:k + 1
        [fi,Dfi] = Legendre(i - 1,-h/2,h/2);
        fj = Legendre(j - 1,-h/2,h/2);
        fjDfi = @(x) Dfi(x).*fj(x);
        G(i,j) = -quadgk(fjDfi,-h/2,h/2) + fj(h/2)*fi(h/2);
        H(i,j) = fj(h/2)*fi(-h/2);
    end
end

I = speye(N);
v = [2:N,1];
W = I(:,v);
A = kron(I,G) - kron(W,H);

%计算初值
for p = 1:(k + 1)*N
    j = mod(p - 1,k + 1);
    i = (p - 1 - j)/(k + 1) + 1;
    fp = Legendre(j,X(i),X(i + 1));
    g = @(x) fp(x).*f(x);
    F(p) = quadgk(g,X(i),X(i + 1));
end

%将初值装填到C的第一列
C = F';
C1 = reshape(C,[k + 1,N]);
C1 = C1';

%开始求解
while T(end) < tb
    t = CFL*h^((k + 1)/kt);
    if T(end) + t >= tb
        t = tb - T(end);
    end
    T = [T;T(end) + t];
    if kt == 3
        C = [C,RK3(A,C(:,end),t)];
    elseif kt == 4
        C = [C,RK4(A,C(:,end),t)];
    end
end

%利用元胞数组储存每一个时间层上的系数矩阵
E = cell(length(T),1);
for j = 1:length(T)
    R = reshape(C(:,j),[k + 1,N]);
    E{j} = R';
end

C = E;

%这里计算Q在每个单元左端点的值
for j = 1:length(T)
    S = C{j};
    Le = cell(k + 1);
    for i = 1:k + 1
        Le{i} = Legendre(i - 1,-h/2,h/2);
    end
    P = zeros(1,N);
    for i = 1:N
        for q = 1:k + 1
            lrd = Le{q};
            P(i) = P(i) + S(i,q)*lrd(-h/2);
        end
    end
    Q(:,j) = P;
end

Q = [Q;Q(1,:)];

end