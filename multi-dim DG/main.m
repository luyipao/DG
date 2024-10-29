% generate x-direction mesh
% some physical constants
n_multiple = 4;
alpha = 2.43694;
c_p = 1830349;
epsilon_r = 11.7;
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
rho_h = @(x) pi * w_max * 2 * N_D(x);

%% step 2 covert equation 3.19 into matrix computation
% P[i, j , 1:2] : I_i, jth basis function value at left and right.
delta_X = X(2:end) - X(1:end-1);
P = zeros(N_x, l + 1, 2);  
P(:, 1, 2) = sqrt(3) ./ sqrt(delta_X');  
P(:, 1, 1) = P(:, 1, 2) .* (-1).^(1:N_x)';
% PP = zeros(N_x, l + 1, 2, N_x, l + 1, 2);
% for i = 1:N_x
%     for j = 1:2
%         for k = 1:2
%             PP(i-1:i+1, j, 1, i, k, 1) = P(i-1:i+1, j, 1) * P(i, k, 1);
%             PP(i-1:i+1, j, 2, i, k , 1) = P(i-1:i+1, j, 2) * P(i, k, 1);
%             PP(i-1:i+1, j, 1, i, k, 2) = P(i-1:i+1, j, 1) * P(i, k, 2);
%             PP(i-1:i+1, j, 2, i, k , 2) = P(i-1:i+1, j, 2) * P(i, k, 2);
%         end
%     end
% end

% generate
A = zeros(2,N_x, 2 * N_x);
B = zeros(2, N_x, 2 * N_x);
C = zeros(2, N_x, 2 * N_x);
D = zeros(2, N_x, 2 * N_x);

for k = 1:2
    for i = 1:N_x-1
        A(k, i,2 * i - 1 : 2 * i + 2) = [-(k == 1) - P(i,1,1) * P(i, k, 2), -(k == 2) - P(i,2,1) * P(i, k, 2), P(i + 1, 1, 1) * P(i, k, 2), P(i + 1, 2, 1) * P(i, k, 2)];
        if i == 1
            B(k, i, 1:4) = [P(i, 1, 2) * P(i, k, 2) + P(i, 1, 1) * P(i, k , 1), P(i, 2, 2) * P(i, k , 2) + P(i, 2, 2) * P(i, k, 2), -P(i + 1, 1, 1) * P(i, k, 2), -P(i + 1, 2, 1)*P(i,k,2)];
            B(k, i, end - 1 : end) = [-P(end, 1,2) * P(i, k, 1), -P(end, 2, 2) * P(i, k, 1)];
        else
            B(k, i, 2*i-3:2*i+2) = [-P(i-1, 1,2) * P(i, k, 1), -P(i-1, 2, 2) * P(i, k, 1),P(i, 1, 2) * P(i, k, 2) + P(i, 1, 1) * P(i, k , 1), P(i, 2, 2) * P(i, k , 2) + P(i, 2, 2) * P(i, k, 2), -P(i + 1, 1, 1) * P(i, k, 2), -P(i + 1, 2, 1)*P(i,k,2)];
        end
        C(k, i, 2 * i - 1: 2 * i) = [(1 == k), 2==k];
        if i == 1
            D(k, 1, 1:2) = [-P(i,1,2)*P(i,k,2), 2*sqrt(3)*(k==2)-P(i,2,2)*P(i,k,2)];
            D(k, 1, end-1:end) =  [P(end,1,2)*P(i,k,1), P(end,2,2)*P(i,k,1)];
        else
            D(k, i, 2*i-3:2*i) = [P(i-1,1,2)*P(i,k,1), P(i-1,2,2)*P(i,k,1), -P(i,1,2)*P(i,k,2), 2*sqrt(3)*(k==2)-P(i,2,2)*P(i,k,2)];
        end
    end
    C(k, N_x, 2 * N_x - 1: 2 * N_x) = [(1 == k), 2==k];
    D(k, end, end-3:end) = [P(end-1,1,2)*P(i,k,1), P(end-1,2,2)*P(i,k,1), 0, 2*sqrt(3)*(k==2)-P(i,2,2)*P(i,k,2)];
    D(k, end, 1:2) = [-P(1,1,1)*P(end,k,2), -P(end,2,1)*P(end,k,2)];
    B(k, end, end - 3:end) = [-P(i-1, 1,2) * P(i, k, 1), -P(i-1, 2, 2) * P(i, k, 1),P(i, 1, 2) * P(i, k, 2) + P(i, 1, 1) * P(i, k , 1), P(i, 2, 2) * P(i, k , 2) + P(i, 2, 2) * P(i, k, 2)];
    B(k, end, 1:2) = [-P(1, 1, 1) * P(i, k, 2), -P(1, 2, 1)*P(i,k,2)];
    A(k, end, end - 1 : end) = [-(k == 1) - P(i, 1, 1) * P(i, k, 1)  + P(i, 1, 2) * P(i, k, 2), -(k == 2) - P(i, 2, 1) * P(i, k, 1) + P(i, 2, 2) * P(i, k, 2)];
end
A = [squeeze(A(1,:,:)), squeeze(B(1,:,:)); 
     squeeze(C(1,:,:)), squeeze(D(1,:,:)); 
     squeeze(A(2,:,:)), squeeze(B(2,:,:)); 
     squeeze(C(2,:,:)), squeeze(D(2,:,:))];


B_str = "c_p(rho_h(x)-N_D(x))";