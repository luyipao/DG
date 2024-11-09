% generate x-direction mesh
% some physical constants
n_multiple = 4;
alpha = 2.43694;
c_p = 1830349;
epsilon_r = 11.7;
c_v = 10;
l = 1; % max dim of basis function
V_bias = 1;
hbar = 1.05457180013e-34;
k_B = 1.380649e-23;
T_L = 300;
m_0 = 9.109383713928e-31;
m_star = 0.32 * m_0;
q = 1.602176634 * 1e-19;
l_star = 1e-6;
t_star = 1e-12;
V_star = 1;
E_star = 0.1  * V_star / l_star;
epsilon_0 = 8.85419e-12;
c_p = (sqrt(2*m_star*k_B*T_L) / hbar)^3 * l_star^2 * q / epsilon_0
c_x = t_star / l_star * sqrt(2 * k_B * T_L / m_star)
c_k = t_star * q * E_star / sqrt(2*m_star * k_B * T_L)
% mesh region
L = 2;



% grid width
delta_x_pos = 0.01;
delta_x_neg = 0.005;
delta_w = alpha / n_multiple;
delta_mu_less0dot7 = 1.7 / 12;
delta_mu_greater0dot7 = 0.3 / 12;

% 
w_max = N_w * delta_w;

% generate grid
X = [0:delta_x_pos:0.2, 0.2:delta_x_neg:0.4, 0.4:delta_x_pos:L];
X = unique(X);
W = linspace(0, w_max, N_w + 1);
MU = [linspace(-1, 0.7, N_mu / 2 + 1), linspace(0.7, 1, N_mu / 2 + 1)];
MU = unique(MU);

% cells num in diff directions
N_x = length(X) - 1;
N_w = 60;
N_mu = 24;

%% step 2 covert equation 3.19 into matrix computation
% P[i, j , 1:2] : I_i, jth basis function value at left and right.
delta_X = X(2:end) - X(1:end-1);
P = zeros(N_x, l + 1, 2);
P(:, 1, 2) = sqrt(1) ./ sqrt(delta_X(:));  
P(:, 1, 1) = P(:, 1, 2) .* (-1).^(0);
P(:, 2, 2) = sqrt(3) ./ sqrt(delta_X(:));
P(:, 2, 1) = P(:, 2, 2) .* (-1).^(1);

% A1 true
A1 = cell(1, N_x);
for i = 1:N_x
    A1{i} = P(i,:,1)' * P(i,:,1);
end
A1 = blkdiag(A1{:});
A1 = A1 + kron(diag(sqrt(3) ./ delta_X), [0 2; 0 0]);
A1 = - epsilon_r * A1;

% dv
A2 = cell(1, N_x);
for i = 1:N_x-1
    A2{i} = epsilon_r * P(i+1,:,1)' * P(i,:,2);
end
A2{N_x} = epsilon_r * P(N_x,:,2)' * P(N_x,:,2);
temp = diag(diag(eye(N_x-1)), 1);
temp(N_x, N_x) = 1;
temp = kron(temp,eye(l+1));
A2 = blkdiag(A2{:});
A2 = A2 * temp;

A = A1 + A2;

B1 = cell(1,N_x);
for i = 1:N_x-1
    B1{i} = -P(i+1,:,1)' * P(i,:,2);
end
B1{N_x} = eye(l+1);
B1 = blkdiag(B1{:});
temp = diag(diag(eye(N_x-1)), 1);
temp = kron(temp,eye(l+1));
B1 = B1 * temp;
b1 = zeros(2*N_x, 1);
b1(end-l:end, 1) = - V_bias * P(end,:,2)';

B2 = cell(1,N_x);
for i = 1:N_x
    B2{i} = P(i,:,2)' * P(i,:,2) + P(i,:,1)' * P(i,:,1);
end
B2 = blkdiag(B2{:});


B3 = cell(1, N_x);
for i = 1:N_x-1
    B3{i}= -P(i,:,2)' * P(i+1,:,1);
end 
B3{N_x} = eye(2);
temp = diag(diag(eye(N_x-1)), -1);
temp = kron(temp, eye(l+1));
B3 = blkdiag(B3{:});
B3 = temp * B3;
b3 = zeros(2*N_x,1);

B = B1 + B2 + B3;

C = eye(2*N_x);

temp = diag(sqrt(3) ./ delta_X);
D1 = kron(temp, [0 2; 0 0]);

D2 = cell(1, N_x);
for i = 1:N_x-1
    D2{i} = -P(i,:,2)' * P(i,:,2);
end
D2{N_x} = -P(1,:,1)' * P(N_x,:,2);
D2 = blkdiag(D2{:});
temp = eye(N_x);
temp(N_x, N_x) = 0;
temp = kron(temp, eye(l+1));
D2 = D2 * temp;
d2 = zeros(2*N_x, 1);
d2(end-l:end,1) = - V_bias * P(end,:,2)';


D3 = cell(1, N_x);
for i = 1:N_x-1
    D3{i} = P(i,:,2)' * P(i+1,:,1);
end
D3{N_x} = eye(2);
D3 = blkdiag(D3{:});
temp = diag(ones(1,N_x-1), -1);
temp = kron(temp, eye(l+1));
D3 = temp * D3;
d3 = zeros(2*N_x, 1);


D = D1 + D2 + D3;


% basis functions
basisFuncs = {@(x) 1 + 0*x, @(x) x};
scr_N_D =@(x) (sqrt(2*0.32*m_0*k_B*T_L) / hbar)^-3 * N_D(l_star * x);

rho_h = @(x) scr_N_D(x);
F1 = arrayfun(@(a, b) gaussLegendre(@(x) c_p * (rho_h(x) - scr_N_D(x)) .* basisFuncs{1}(((b-a)*x+a+b)/2), a, b), X(1:end-1), X(2:end));
F1 = F1 ./ sqrt(delta_X);
F2 = arrayfun(@(a, b) gaussLegendre(@(x) c_p * (rho_h(x) - scr_N_D(x)) .* basisFuncs{2}(((b-a)*x+a+b)/2), a, b), X(1:end-1), X(2:end));
F2 = F2 * sqrt(3) ./ sqrt(delta_X);

F = zeros(length(F1) + length(F2), 1);
F(1:2:end) = F1;
F(2:2:end) = F2;
F1 = F - b1 - b3;
F2 = - d2 - d3;


temp = [A B; C D];
coeff = [A B; C D] \ [F1; F2];
C_q = coeff(1:2*N_x,1);
C_Psi = coeff(2*N_x + 1:end,1);
Psi_h = basisPolys(X, reshape(coeff(2*N_x+1:end), l + 1, []), l, basisFuncs, [0, 1]');
E = basisPolys(X, reshape(coeff(1:2*N_x), l + 1, []), l, basisFuncs, [0, 1]');
% draw p_h and Psi_h
XX = linspace(0, L,10000);
clf
hold on
plot(XX, Psi_h.solve(XX));
scatter(XX, E.solve(XX), '.');
legend("Psi_h", "E");
hold off