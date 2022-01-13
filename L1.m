function res = L1(v)
% L1 Summary of this function goes here
%   Detailed explanation goes here
    v1 = v(1:3, 1);
    v2 = v(4:6, 1);
    res = zeros(6, 6);
    res(1:3, 1:3) = hat(v2);
    res(4:6, 1:3) = hat(v1);
    res(4:6, 4:6) = hat(v2);
end