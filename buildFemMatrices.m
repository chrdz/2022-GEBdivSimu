function [M, K, P1, P2, P3, Pdagger1, Pdagger2, Pdagger3, LL0, LL1] = buildFemMatrices(JJ, AA, BB, GG, GGdagger, Q_lyap, NNB, he, Ne, massMat, flexMat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% JJ y_t + AA y_x + BB y = GG(y)y
% G is d_G (diagonal system) or G (physical system)
    
    Ntot = numel(NNB);
    Ni = 12;

    %% Element matrices
    syms XI;
    Ntild =  [(1-XI)*(1-2*XI), 4*XI*(1-XI), XI*(2*XI - 1)];
    Dxi_Ntild = diff(Ntild, XI);

    Me_syms = int((Ntild')*Ntild, 0, 1);
    Ke_syms = int((Dxi_Ntild')*Ntild, 0, 1);
    P1e_syms = int(Ntild(1)*(Ntild')*Ntild, 0, 1);
    P2e_syms = int(Ntild(2)*(Ntild')*Ntild, 0, 1);
    P3e_syms = int(Ntild(3)*(Ntild')*Ntild, 0, 1);
    LL1e_syms = int((Dxi_Ntild')*Dxi_Ntild, 0, 1);

    Me = double(Me_syms);
    Ke = double(Ke_syms);
    P1e = double(P1e_syms);
    P2e = double(P2e_syms);
    P3e = double(P3e_syms);
    LL1e = double(LL1e_syms);

    %% Assemble the mass and stiffness matrices
    % and matrices LL0 and LL1 for the lyap func
    disp('Assembling the mass and stiffness matrices..')
    tic

    M = sparse(Ntot,Ntot);    % initialize zero mass matrix 
    K = sparse(Ntot,Ntot);    % initialize zero stiffness matrix
    LL0= sparse(Ntot,Ntot);
    LL1= sparse(Ntot,Ntot);
    

    % integrals
    for ii = 1:Ni
        for jj = 1:Ni
            for ee = 1:Ne    
                idxR = [NNB(ii, 2*ee-1), NNB(ii, 2*ee), NNB(ii, 2*ee+1)];
                idxC = [NNB(jj, 2*ee-1), NNB(jj, 2*ee), NNB(jj, 2*ee+1)];
                M(idxR, idxC) = M(idxR, idxC) + JJ(ii, jj)*he*Me;
                K(idxR, idxC) = K(idxR, idxC) + BB(ii, jj)*he*Me - AA(ii, jj)*Ke;
                LL0(idxR, idxC) = LL0(idxR, idxC) + Q_lyap(ii, jj)*he*Me;
                LL1(idxR, idxC) = LL1(idxR, idxC) + Q_lyap(ii, jj)/he*LL1e;
            end
        end
    end

    toc

    %% Assemble the matrices that define the nonlinearity
    disp('Assembling the matrices defining the map Q(y)..')
    tic

    P1 = cell(Ni, Ne);
    P2 = cell(Ni, Ne);
    P3 = cell(Ni, Ne);
    Pdagger1 = cell(Ni, Ne);
    Pdagger2 = cell(Ni, Ne);
    Pdagger3 = cell(Ni, Ne);
    for pp = 1:Ni
        ep = zeros(12, 1); ep(pp, 1) = 1;
        Gp = GG(ep, massMat, flexMat);
        Gdaggerp = GGdagger(ep, massMat, flexMat);
        for ee = 1:Ne
            P1pe = sparse(Ntot, Ntot);   % Initialize zero matrices
            P2pe = sparse(Ntot, Ntot);
            P3pe = sparse(Ntot, Ntot);
            Pdagger1pe = sparse(Ntot, Ntot);
            Pdagger2pe = sparse(Ntot, Ntot);
            Pdagger3pe = sparse(Ntot, Ntot);
            for ii=1:Ni
                for jj = 1:Ni
                    idxR = [NNB(ii, 2*ee-1), NNB(ii, 2*ee), NNB(ii, 2*ee+1)];
                    idxC = [NNB(jj, 2*ee-1), NNB(jj, 2*ee), NNB(jj, 2*ee+1)];
                    P1pe(idxR, idxC) = P1pe(idxR, idxC) + Gp(ii, jj)*he*P1e;
                    P2pe(idxR, idxC) = P2pe(idxR, idxC) + Gp(ii, jj)*he*P2e;
                    P3pe(idxR, idxC) = P3pe(idxR, idxC) + Gp(ii, jj)*he*P3e;
                    Pdagger1pe(idxR, idxC) = Pdagger1pe(idxR, idxC) + Gdaggerp(ii, jj)*he*P1e;
                    Pdagger2pe(idxR, idxC) = Pdagger2pe(idxR, idxC) + Gdaggerp(ii, jj)*he*P2e;
                    Pdagger3pe(idxR, idxC) = Pdagger3pe(idxR, idxC) + Gdaggerp(ii, jj)*he*P3e;
                end
            end
            P1(pp, ee) = {P1pe};
            P2(pp, ee) = {P2pe};
            P3(pp, ee) = {P3pe};  
            Pdagger1(pp, ee) = {Pdagger1pe};
            Pdagger2(pp, ee) = {Pdagger2pe};
            Pdagger3(pp, ee) = {Pdagger3pe}; 
        end
    end

    toc


end

