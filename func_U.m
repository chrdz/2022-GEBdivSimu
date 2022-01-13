function res = func_U(v)
% linear function used to transform the pde with unknown rotation matrix R
% to a PDE with unknown the quaternion that parametrizes R
    res = zeros(4, 4);
    res(1, 2:4) = -transpose(v)/2;
    res(2:4, 1) = v/2;
    res(2:4, 2:4) = -hat(v)/2;
end