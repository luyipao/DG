%% 
f = @(t) 2 * pi * cos(2* pi .*t);


%%  main 
u = @(t) sin(2*pi .* t);
a = 0;
b = 2;
initial_value = 0;
for p = [1,inf]
    for degree = [0,1,2]
        for size = [0.1, 1/20, 1/40, 1/80,1/160, 1/ 320]
            [~,C] = DG(f, a,b,initial_value, size,degree);
            u_h = general(C,a,b,size);
            e_h = error_vec(u,u_h,a,b);
            disp(['h is ', num2str(size), '; k is: ', num2str(degree), '; in integral [', num2str(a), ', ', num2str(b), ' ], norm is : ', num2str(p), ', e_h is ', num2str(h * norm(e_h,p))]);
        end
    end
end

function u_hh = general(C,xa,xb,h)
X = xa:h:xb;
N = length(X) - 1;
u_h = cell(1, N);
for i = 1:N
    u_h{i} = @(t) 0;
end
u_hh = @(t) 0;
[K N] = size(C);
for ii = 1:N
    for jj = 1:K
        [phi,~] = base_function(jj-1,X(ii),X(ii+1));
        u_h{ii} = @(t) u_h{ii}(t) + C(jj,ii) * phi(t);
    end
    u_hh  = @(t) u_hh(t) + (X(ii) < t) .*(t <= X(ii+1)) .* u_h{ii}(t);
end
u_hh = @(t) u_hh(t) + (t == xa) .* u_h{1}(t);
end

function [phi, diff_phi] = base_function(k,a,b)
size = b - a;
phi = @(x) ((x-a).^k / (size^k)) * (k>=1) + 1.0 * (k ==0);
% derivation of base function psizei
diff_phi = @(x) (k >=1) * (k * (x-a).^(k-1) / size^k);
end

%% Function
%figure
% input: u真实函数；u_h数值解函数，[a,b]网格区间
% output: 真实函数图像和数值解函数图像
 function figure(u,u_h,a,b)
    x = a:0.001:b;
    plot(x,u(x),'g',x,u_h(x),'--');
    legend('exact solution','numericial solution');
 end
 % e_h
 % input: u真实函数；u_h数值解函数，[a,b]网格区间
 % output: e_h误差向量
 function e_h = error_vec(u,u_h,a,b)
    x = a:0.001:b;
    e_h = u(x) - u_h(x);
 end