%%
f = @(t) cos(t);
% base function psizei
base_func = @(t,xi,j,size) ((t-xi).^j / (size.^j)) .* (j>=1) + 1.0 .* (j ==0) ;
% derivation of base function psizei
diff_base_func = @(t,xi,j,size) (j >=1) .* (j .* (t-xi).^(j-1) / size^j) ;


%% 
u = @(t) sin(t);
a = 0;
b = 2;
initial_value = 0;
for p = [1,inf]
    for degree = [0,1,2]
        for size = [0.1, 1/20, 1/40, 1/80,1/160, 1/ 320]
            u_h = DG(f,  a,b,initial_value, size,degree,base_func,diff_base_func);
            e_h = error_vec(u,u_h,a,b);
            disp(['h is ', num2str(size), '; k is: ', num2str(degree), '; in integral [', num2str(a), ', ', num2str(b), ' ], norm is : ', num2str(p), ', e_h is ', num2str(norm(e_h,p))]);
        end
    end
end

%% Function
% input: size间隔；degree基函数的阶；区间[a,b]；f；initial_value初值;base_func
% 基函数；diff_base_func基函数的导数
% output: u_h数值解
% tdl: 该DG于base绑定，只是适用当前base function
function u_hh = DG(f,lEndpoint,rEndpoint,initial_value,size,degree,base_func,diff_base_func) 
N = fix((rEndpoint-lEndpoint)/size);
A = zeros(degree+1,degree+1);
b = zeros(degree+1,1);
u_h = cell(1, N);
x = lEndpoint:size:rEndpoint;
u_hh = @(t) 0;
for i = 1:N
    u_h{i} = @(t) 0;
end

for ii = 1:N
for m = 1:degree+1
    for j = 1:degree+1
        A(m,j) = -1.0 .* integral(@(t) base_func(t,x(ii),j-1,size).*diff_base_func(t,x(ii), m-1, size),  x(ii),  x(ii+1)) + base_func(x(ii+1),x(ii),j-1,size).*base_func(x(ii+1),x(ii),m-1,size);
    end
    b(m) = initial_value .* base_func(x(ii),x(ii),m-1,size) + integral(@(t) base_func(t,x(ii),m-1,size).*f(t),x(ii),x(ii+1));
end
c = A \ b;
for ll = 1:(degree+1)
    u_h{ii} = @(t) u_h{ii}(t) +  c(ll) .* base_func(t,x(ii),ll-1,size);
end
if (ii ~= N)
    initial_value = u_h{ii}(x(ii+1));
end
end
for jj = 1:N
    u_hh = @(t) u_hh(t) + ( x(jj) <= t) .*(t < x(jj+1)) .* u_h{jj}(t);
end
u_hh = @(t) u_hh(t) + (t == rEndpoint) .* u_h{N}(t);
end
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