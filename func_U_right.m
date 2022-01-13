function res = func_U_right(v)
% linear function similar to U, used in the paper of Zupan
% we probably will remove this function if we can't make work the code of Zupan
    res = zeros(4, 4);
    res(1, 2:4) = -transpose(v);
    res(2:4, 1) = v;
    res(2:4, 2:4) = hat(v);
end