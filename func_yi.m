function res = func_yi(x, i, n, Y, Nx, NNB,Ne,he, xList)
% func_y return the real FEM approx of y
% x in [0, \ell] is the spatial variable
% n is the time index
% NNB is in the numbering of the FEM system
    res = 0;
    for kk = 1:Nx
        res = res + Nk(x, kk,Ne,he, xList)*Y(NNB(i, kk), n);
    end
end