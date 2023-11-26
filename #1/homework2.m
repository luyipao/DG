f = @(x) 2 * pi * cos(2* pi .*x);
% u(x,t) = sin(x - t)
%input: 区间[xa,xb]；cell宽度size；解空间多项式最高阶degree；t=0时函数u(x,0)=u_0
% output: 系数矩阵C，其列向量c_i表示在cell I_i上的u的系数
u_0 =@(t)  sin(2*pi .* t);
u = @(x,t) sin(2*pi .*(x-t));
for h = [1/10,1/20,1/40]
     degree = 2;
    [u_hh,U_h] = RKDG(f,0,1,0,1,u_0,h,degree);
    x = 0:h:1;
    t = 0:h:1;
    U = zeros(1/h+1,1/h+1);
    for ii = 1:1/h+1
        for tt = 1:1/h+1
            U(tt,ii) = u(x(ii),t(tt));
        end
    end
    E_h = U - U_h;
    disp(['h is ', num2str(h), '; k is: ', num2str(degree), '; error is ', num2str(norm(E_h,1))]);
end
function [u_hh,D] = RKDG(f,xa,xb,ta,tb,u_0,size,degree)
    X = xa:size:xb;
    N = length(X) - 1;
    h_t = (tb-ta) / length(X);
    A = zeros(degree + 1, degree+1);
    u_hh = cell(length(X),1);% t_n 上近似[xa,xb]的u_h
    [u_h,C] = DG(f,xa,xb,u_0(xa),size,degree);
    for tt = 1:N
        u_hh{tt} = builder(C,xa,xb,size);
        for ii = 1:N
            for m = 1:degree+1
                for jj = 1:degree+1  
                    [phi_m,~] = base_func(m,X(ii),X(ii+1));
                    [phi_jj,~] = base_func(jj,X(ii),X(ii+1));
                    A(m,jj) = integral(@(x) phi_m(x) .* phi_jj(x), X(ii), X(ii+1));
                end
            b = L(u_h, X(ii), X(ii+1), degree + 1);    
            end
            C(:,ii) = C(:,ii) + h_t  * (A\b);
        end
    end
    u_hh{length(X)} = builder(C,xa,xb,size);
    D = zeros(length(X),length(X));
    for tt = 1:length(X)
        D(tt,:) = u_hh{tt}(X);
    end
end


%L base on base_function, 
% input: t_n时初值u^{(n)}，时间间隔t，算子L(x)
% ouput: u^{(n+1)}

function u = RK3(u,t,L,xa,xb,K)
u1 = u + t * L(u,xa,xb,K);
u2 = 3 / 4 * u + 1/4 * (u1 + t .* L(u1));
u = 1/3 * u + 2/3 * (u2 + t .* L(u2));
end


% input: 区间[a,b]，阶k
% output: phi [a,b]上基函数; diff_phi phi的导函数
function u_hh = builder(C,xa,xb,h)
X = xa:h:xb;
N = length(X) - 1;
u_h = cell(1, N);
for i = 1:N
    u_h{i} = @(t) 0;
end
u_hh = @(t) 0;
[K, ~] = size(C);
for ii = 1:N
    for jj = 1:K
        [phi,~] = base_func(jj-1,X(ii),X(ii+1));
        u_h{ii} = @(t) u_h{ii}(t) + C(jj,ii) * phi(t);
    end
    u_hh  = @(t) u_hh(t) + (X(ii) < t) .*(t <= X(ii+1)) .* u_h{ii}(t);
end
u_hh = @(t) u_hh(t) + (t == xa) .* u_h{1}(t);
end

% input: 阶k, cell [a,b]
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
    [phi,diff_phi] = base_func(ii,xa,xb);
    b(ii) = integral(@(x) u(x) .* diff_phi(x), xa, xb) - u(xb) * phi(xb) + u(xa) * phi(xa);
end
end