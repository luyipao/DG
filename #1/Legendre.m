% Function: Legendre
% Input: Degree 'n' as a vector and interval [a,b]
% Input condition: b > a
% Output: Two sym vectors 'PN' & 'DPN', where 'DPN' is the derivative of 'PN'
% Output condition: Each cell contains the corresponding base function


function [PN, DPN] = Legendre(n, a, b)
    syms x; % defining the symbolic variable

    % Constructing Legendre Polynomials
    P = legendreP(n, (2*x-a-b)/(b-a));
    % Getting the coefficients
    c = coeffs(P);
    % Making the polynomials monic
    PN = P/c(end); % use parentheses instead of curly braces
    % Calculating the derivative
    DPN = diff(PN); % use parentheses instead of curly braces
    % sym to function handle
    PN = matlabFunction(PN, 'Vars', x);
    DPN = matlabFunction(DPN, 'Vars', x);
end


% Example usage - considering the list of degrees [0, 1, 2]
% [PN, DPN] = Legendre([0,1,2], -1, 1);
% [PN, DPN] = Legendre(2,-1,1);

