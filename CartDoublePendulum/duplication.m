function [D] = duplication(n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% n: the size of the full square matrix. i.e. if 
% the matrix is a 2x2, then n = 2;


m = n * (n + 1)/2;  % Number of elements in the half-vectorization
nsq = n ^ 2;        % Total number of elements in the matrix

r = 1;
a = 1;
v = zeros(1,nsq);

% Create an index set for the sparse duplication matrix
for i = 1:n
    v(r:r + i - 2) = i - n + cumsum(n - (0:i-2));
    r = r + i - 1;
    
    v(r:r + n - i) = a:a + n - i;
    r = r + n - i + 1;
    a = a + n - i + 1;    
end
% Create the duplication matrix as a sparse matrix
D = sparse(1:nsq, v, 1, nsq, m);

end

