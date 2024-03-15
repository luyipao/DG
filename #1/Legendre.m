% Function: Legendre
% Input: Degree 'n' as a vector and interval [a,b]
% Input condition: b > a
% Output: Two sym vectors 'PN' & 'DPN', where 'DPN' is the derivative of 'PN'
% Output condition: Each cell contains the corresponding base function


function [PN, DPN] = Legendre(n, a, b)
    syms x; % defining the symbolic variable
    nn = numel(n);
    PN = sym('PN', [nn 1]); % Now PN and DPN are symbolic vectors
    DPN = sym('DPN', [nn 1]);

    for i = 1:numel(n) % looping over every element in 'n'
        % Constructing Legendre Polynomials
        P = legendreP(n(i), (2*x-a-b)/(b-a));
        % Getting the coefficients
        c = coeffs(P);
        % Making the polynomials monic
        PN(i) = P/c(end); % use parentheses instead of curly braces
        % Calculating the derivative
        DPN(i) = diff(PN(i)); % use parentheses instead of curly braces
    end
    % sym convert to function
    PN = matlabFunction(PN);
    DPN = matlabFunction(DPN);
end


% Example usage - considering the list of degrees [0, 1, 2]
% [PN, DPN] = Legendre([0,1,2], -1, 1);
% [PN, DPN] = Legendre(2,-1,1);