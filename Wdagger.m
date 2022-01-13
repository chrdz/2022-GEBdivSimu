function res = Wdagger(i, v)
% linear function used to compute the Jacobian of W
    
if i == 0
    res = 2*[v, hat(v)];
else
    ei = zeros(3, 1); ei(i, 1) = 1;
    res = 2*[hat(v)*ei, ei*transpose(v) - v*transpose(ei) + v(i, 1)*eye(3)];
end