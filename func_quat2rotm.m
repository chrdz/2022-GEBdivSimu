function res = func_quat2rotm(quat)
% transfor a unit norm quaternion to the corresponding rotation matrix
    q0 = quat(1, 1);
    q = quat(2:4, 1);
    res = (q0^2 - transpose(q)*q)*eye(3) + 2*q*transpose(q) + 2*q0*hat(q);
end