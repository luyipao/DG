% generate x-direction mesh
% some physical constants
n_multiple = 4;
alpha = 2.43694;
c_p = 1830349;
epsilon_r = 11.7;
c_v = 10;
l = 1; % max dim of basis function
% mesh region
L = 1;

% cells num in diff directions
N_x = 120;
N_w = 60;
N_mu = 24;

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

% setp 1
rho_h = @(x) pi * w_max * 2 * N_D(x)  ;

%% step 2 covert equation 3.19 into matrix computation
% P[i, j , 1:2] : I_i, jth basis function value at left and right.
delta_X = X(2:end) - X(1:end-1);
P = zeros(N_x, l + 1, 2);
P(:, 1, 2) = sqrt(1) ./ sqrt(delta_X(:));  
P(:, 1, 1) = P(:, 1, 2) .* (-1).^(1:N_x)';
P(:, 2, 2) = sqrt(3) ./ sqrt(delta_X(:));
P(:, 2, 1) = P(:, 2, 2) .* (-1).^(1:N_x)';

A1 = cell(1, N_x);
for i = 1:N_x
    A1{i} = P(i,:,1)' * P(i,:,1);
end
A1 = blkdiag(A1{:});
A1 = A1 + kron(diag(sqrt(3) ./ delta_X), [0 2; 0 0]);
A1 = -epsilon_r * A1;

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
B1{N_x} = -P(1,:,1)' * P(N_x,:,2);
B1 = blkdiag(B1{:});
temp = diag(diag(eye(N_x-1)), 1);
temp(N_x, 1) = 1;
temp = kron(temp,eye(l+1));
B1 = B1 * temp;
clear temp;

B2 = cell(1,N_x);
for i = 1:N_x
    B2{i} = P(i,:,2)' * P(i,:,2) + P(i,:,1)' * P(i,:,1);
end
B2 = blkdiag(B2{:});
B2 = kron(eye(N_x), eye(l+1)) * B2;

B3 = cell(1, N_x);
for i = 2:N_x
    B3{i}= -P(i-1,:,2)' * P(i,:,1);
end
B3{1} = -P(N_x,:,2)' * P(1,:,1);
temp = diag(diag(eye(N_x-1)), -1);
temp(1,N_x) = 1;
temp = kron(temp, eye(l+1));
B3 = blkdiag(B3{:});
B3 = temp * B3;

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
temp(N_x, 1) = 1;
temp = kron(temp, eye(l+1));
D2 = D2 * temp;

D3 = cell(1, N_x);
for i = 2:N_x
    D3{i} = P(i-1,:,2)' * P(i,:,1);
end
D3{1} = P(N_x,:,2)' * P(1,:,1);
D3 = blkdiag(D3{:});
temp = diag(ones(1,N_x-1), -1);
temp(1,N_x) = 1;
temp = kron(temp, eye(l+1));
D3 = temp * D3;

D = D1 + D2 + D3;


% basis functions
basisFuncs = {@(x) 1 + 0*x, @(x) x};

F1 = arrayfun(@(a, b) gaussLegendre(@(x) c_p * (rho_h(x) - N_D(x)) .* basisFuncs{1}(((b-a)*x+a+b)/2), a, b), X(1:end-1), X(2:end));
F1 = F1 ./ sqrt(delta_X);
F2 = arrayfun(@(a, b) gaussLegendre(@(x) c_p * (rho_h(x) - N_D(x)) .* basisFuncs{2}(((b-a)*x+a+b)/2), a, b), X(1:end-1), X(2:end));
F2 = F2 * sqrt(3) ./ sqrt(delta_X);

F = zeros(length(F1) + length(F2), 1); 
F(1:2:end) = F1; 
F(2:2:end) = F2; 
F = [F(:); zeros(2*N_x, 1)];

coeff = [A B; C D] \ F;








p_h = basisPolys(X, reshape(coeff(1:2*N_x), l + 1, []), l, basisFuncs, [0, 1]');
Psi_h = basisPolys(X, reshape(coeff(2*N_x+1:end), l + 1, []), l, basisFuncs, [0, 1]');
E = basisPolys(X, reshape(-c_v * coeff(1:2*N_x), l + 1, []), l, basisFuncs, [0, 1]');
% draw p_h and Psi_h
X = linspace(0, 1,10000);
clf
hold on
plot(X, p_h.solve(X));
plot(X, Psi_h.solve(X));
plot(X, E.solve(X));
legend("p_h", "Psi_h", "B2");
hold off