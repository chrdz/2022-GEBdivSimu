function res = hat(u)
% hat Summary of this function goes here
%   Detailed explanation goes here
    res = zeros(3, 3);
    res(1, 2) = -u(3, 1); res(1, 3) = u(2, 1);
    res(2, 1) = u(3, 1); res(2, 3) = -u(1, 1);
    res(3, 1) = -u(2, 1); res(3, 2) = u(1, 1);
end