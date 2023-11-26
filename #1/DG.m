% input: size cell长度；degree基函数的阶；区间[a,b]；f；initial_value初值
% output: u_hh数值解；C系数矩阵，C(:,i)表示u在I_i上关于基函数的系数
% tdl: 该DG于base绑定，只是适用当前base function
function [u_hh,C] = DG(f,lEndpoint,rEndpoint,initial_value,size,degree)
base_func = @(t,xi,j,size) ((t-xi).^j / (size.^j)) .* (j>=1) + 1.0 .* (j ==0);
% derivation of base function psizei
diff_base_func = @(t,xi,j,size) (j >=1) .* (j .* (t-xi).^(j-1) / size^j);
N = fix((rEndpoint-lEndpoint)/size);
A = zeros(degree+1,degree+1);
b = zeros(degree+1,1);
u_h = cell(1, N);
for i = 1:N
    u_h{i} = @(t) 0;
end
C = zeros(degree + 1, N);
x = lEndpoint:size:rEndpoint;
u_hh = @(t) 0;

for ii = 1:N
for m = 1:degree+1
    for j = 1:degree+1
        A(m,j) = -1.0 .* integral(@(t) base_func(t,x(ii),j-1,size).*diff_base_func(t,x(ii), m-1, size),  x(ii),  x(ii+1)) + base_func(x(ii+1),x(ii),j-1,size).*base_func(x(ii+1),x(ii),m-1,size);
    end
    b(m) = initial_value .* base_func(x(ii),x(ii),m-1,size) + integral(@(t) base_func(t,x(ii),m-1,size).*f(t),x(ii),x(ii+1));
end
c = A \ b;
C(:,ii) = c;
for ll = 1:(degree+1)
     u_h{ii} = @(t) u_h{ii}(t) +  c(ll) .* base_func(t,x(ii),ll-1,size);
 end
if (ii ~= N)
    initial_value = u_h{ii}(x(ii+1));
end
end
% 得到u_hh
for jj = 1:N
    u_hh = @(t) u_hh(t) + ( x(jj) < t) .*(t <= x(jj+1)) .* u_h{jj}(t);
end
u_hh = @(t) u_hh(t) + (t == lEndpoint) .* u_h{1}(t);
end