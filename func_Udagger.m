function res = func_Udagger(qbold)
% linear function used to compute the Jacobian of R
    q0 = qbold(1, 1);
    q = qbold(2:4, 1);
    res = 1/2*[-transpose(q);...
        q0*eye(3)+hat(q)];
end