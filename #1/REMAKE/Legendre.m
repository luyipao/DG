% Function: Legendre
% Input: Degree 'n' as a vector and interval [a,b]
% Input condition: b > a
% Output: Two sym vectors 'PN' & 'DPN', where 'DPN' is the derivative of 'PN'
% Output condition: Each cell contains the corresponding base function


function [PN, DPN] = Legendre(n,a,b)
    P0 = @(x) 1;
    DP0 = @(x) 0;
    P1 = @(x) x;
    DP1 = @(x) 1;
    if n==0
        P = @(x) 1;
        DP = @(x) 0;
    elseif n == 1
        P = @(x) x;
        DP = @(x) 1;
    else
        for ii = 2:n
            P = @(x) ((2.*ii - 1) * x * P1(x) - (ii - 1) * P0(x)) / ii;
            DP = @(x) ((2*ii - 1) * P1(x) + (2*ii - 1) * x * DP1(x) - (ii - 1) * DP0(x)) / ii;
            P0 = P1;
            DP0 = DP1;
            P1 = P;
            DP1 = DP;
        end
    end
    PN = @(x) sqrt((2*n + 1) / (b - a)) * P((2*x - b - a) / (b - a));
    DPN = @(x) sqrt((2*n + 1) / (b - a)) * DP((2*x - b - a) / (b - a));
end



