function res = func_vec(M)
% hat Summary of this function goes here
%   Detailed explanation goes here
    res = zeros(3, 1);
    res(1, 1) = M(3, 2);
    res(2, 1) = M(1, 3);
    res(3, 1) = M(2, 1);
end