%%
% Simulations for a geometrically exact beam 
% the beam is free at one tip (x=ell)
% forces/moments are applied at the other tip (x=0)
% We work with the Intrinsic formulation of the GEB model (Hodges 2003)
% We use P2 elements for the spatial discretization

% We diagonalise the system before solving, 
% after solving we recover the physical unknown and plot it
% The aim is also to have to option of working with the physial system

% Here we additionnaly do a new fixed point, supposedly more precise:
% i.e. during the fixed point to recover y, we recover the orientation of
% the cross sections at x=0 as well.

% Moreover we use the 'p2 method' to recover the position of the beam

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SPACE DISCRETIZATION                            %
% nodes:      1   2   3   4   5     2e-1 2e  2e+1               Nx-1  Nx  %
%             |---o---|---o---|  ...  |---o---|---o---|---o---|---o---|   %
% elements:       1       2               e              Ne-1     Ne      %
%                                                                         %
%                        TIME DISCRETIZATION                              %
% time instances:     1       2                     Nt-1     Nt           %
%                     |-------|-------|  ...  |-------|-------|           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% THE PDE SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We consider the IGEB model.                                            %
% The beam is clamped at x=0.                                            %
% A velocity feedback is applied at x=\ell.                              %
%                                                                        %   
% The unknown is:   y = (v^T, z^T)^T    where:                           %
% v = linear and angular velocities (expressed in body-attached basis),  %
% z = forces and moments (expressed in body-attached basis).             %
%                                                                        %
% The system reads:                                                      %
%                                                                        %
% J(x) y_t + A y_x + B(x)y = g(x,y)        in (0, \ell) x (0, T)         %
% z(0, t) = F(t)                           t in (0, T)                   %
% z(\ell, t) = 0                           t in (0, T)                   %
% y(x,0) = y^0(x)                          x in (0, \ell).               %
%                                                                        %
% or                                                                     %
%                                                                        %
% J(x) y_t + A y_x + B(x)y = g(x,y)        in (0, \ell) x (0, T)         %
% v(0, t) = F(t)                           t in (0, T)                   %
% z(\ell, t) = 0                           t in (0, T)                   %
% y(x,0) = y^0(x)                          x in (0, \ell).               %
%                                                                        %
% M,C are the mass and flexibility matrices.                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all
close all
clc
fs = 14; % font size
size_line = 2; % line size for plots;

%% To be chosen
linearized = false;     % true: we take the IGEB system without g(y)
                        % false: we take the real (semilinear) IGEB system
plot_y = true;     % true: we plot the solution y
                   % false: we don't plot the solution y
plot_centerline = true;   % true: we plot the centerline of the beam
                          % false: we don't plot the centerline
                          
plot_norm = false;

centerline_t_scheme = 0;     % 0: mid point rule
                             % 1: explicit Euler      (does not work properly)
                             % 2: ode45 (RK)
                             % 3: implicit Euler
                             % 4: zupan paper scheme  (does not work properly)           
centerline_x_scheme = 0;

problem = 1;    % 0: flying spaghetti problem 2D
                % 1: flying spaghetti problem 3D
                % 2: book
                
approx_rot = 0;    % 0 most precise mid point everywhere
                   % 1 no mid point for Rfext in W1m
                   % 2 no mid point for Rfext in W1m and in calU(v2)
                
diagonal = true;
                             
ell = 10;    % length of the space interval
Ne = 20;     % number of elements
T = 10;      % end of the time interval
ht = 0.1;   % time step

%% curvature before deformation
kap = [0; 0; 0];  % curvature before deformation
boldE = zeros(6, 6); boldE(5, 3) = -1; boldE(6, 2) = 1;   % initial strain matrix

%% beam parameters
% ---------------------------------------- %
% parameters from the paper:
% Consistent structural linearisation in flexible-body dynamics with large rigid-body motion
% Henrik Hesse, Rafael Palacios, 2012
% ---------------------------------------- %
EA = 10^4; GAs = 10^4;
if problem == 0 | problem == 1
    EI = 500; GJ = 500;
elseif problem == 2
    EI = 1000; GJ = 1000;
end
rhoA = 1; rhoJ = diag([20, 10, 10]);
massMat = blkdiag(rhoA*eye(3), rhoJ);            % the MASS matrix
flexMat = inv(diag([EA, GAs, GAs, GJ, EI, EI])); % the FLEXIBILITY matrix

%% Coeffficients of the PDE system
J = blkdiag(massMat, flexMat);
B = zeros(12, 12); B(1:6, 7:12) = - boldE; B(7:12, 1:6) = transpose(boldE);
A = zeros(12, 12); A(1:6, 7:12) = -eye(6); A(7:12, 1:6) = -eye(6);
Ni = 12;      % number of PDEs

%% Coeffficients of the PDE system DIAGONAL
if diagonal == true
    d_L = 1/sqrt(2)*[[eye(6), eye(6)]; [eye(6), -eye(6)]];
    d_Linv = d_L;
    d_J = d_L*J*d_Linv;
    d_D = [[-eye(6), zeros(6)]; [zeros(6), eye(6)]];
    d_G = @(r, M, C) d_L*G(d_Linv*r, M, C)*d_Linv;
    d_Gdagger = @(r, M, C) d_L*Gdagger(d_Linv*r, M, C)*d_Linv;
    d_B = d_L*B*d_Linv;
    if problem == 0 | problem == 1
        d_W1 = d_D*[[zeros(6), eye(6)]; [zeros(6), eye(6)]];
        d_W2 = -d_D*[[eye(6), zeros(6)]; [eye(6), zeros(6)]];
        d_W3 = sqrt(2)*[zeros(6); eye(6)];
    elseif problem ==  2
        d_W1 = d_D*[[zeros(6), eye(6)]; [zeros(6), eye(6)]];
        d_W2 = -d_D*[[eye(6), zeros(6)]; [-eye(6), zeros(6)]];
        d_W3 = -sqrt(2)*[zeros(6); eye(6)];
    end
end

%% Energy matrix
QP_lyap = blkdiag(massMat, flexMat);
Q_lyap = QP_lyap; % pr moment on regarde juste la norme H1

%% Some dependent variables
Nx = 2*Ne + 1;               % number of nodes
Nt = T/ht +1;                % i.e., ht = T/(Nt-1)
he = ell/Ne;                 % length of one element
x = linspace(0,ell,Nx);      % spatial grid with the node positions
hx = x(2) - x(1);            % spatial step
Ntot = Ni*Nx;                % number of unknowns without BC
t = linspace(0, T, Nt);      % time instances

if linearized == true
    linNonlin = 'lin';
else
    linNonlin = 'nonlin';
end
%file_name_y = ['fig/SOLY_len', num2str(ell), '_', linNonlin, '.pdf'];
%name_ini = ['fig/Y0_len', num2str(ell), '.pdf'];

%% FEM matrices
NNB = reshape(1:Ntot, Ni, Nx);   % node numbers
                                 % i.e. NNB(i, k) = 2*(k-1)+i
                                 
if diagonal ==  true
    [M, K, P1, P2, P3, Pdagger1, Pdagger2, Pdagger3, LL0, LL1] =...
        buildFemMatrices(d_J, d_D, d_B, d_G, d_Gdagger, Q_lyap, ...
        NNB, he, Ne, massMat, flexMat);
else
    [M, K, P1, P2, P3, Pdagger1, Pdagger2, Pdagger3, LL0, LL1] =...
        buildFemMatrices(J, A, B, @G, @Gdagger, Q_lyap, ...
        NNB, he, Ne, massMat, flexMat);
end

%% add boundary terms
U = sparse(Ntot, 6);

if diagonal == true
    for ii = 1:Ni
        for jj = 1:Ni
            K(NNB(ii, Nx), NNB(jj, Nx)) = K(NNB(ii, Nx), NNB(jj, Nx)) + d_W1(ii, jj);
            K(NNB(ii, 1), NNB(jj, 1)) = K(NNB(ii, 1), NNB(jj, 1)) + d_W2(ii, jj);
        end
        for jj = 1:6
            U(NNB(ii, 1), jj) = U(NNB(ii, 1), jj) + d_W3(ii, jj);
        end
    end
else
    if problem == 0 | problem == 1
        for ii = 1:6
            U(NNB(ii, 1), ii) = U(NNB(ii, 1), ii) + 1;
            K(NNB(ii+6, 1), NNB(ii, 1)) = K(NNB(ii+6, 1), NNB(ii, 1)) + 1;
        end
    elseif problem == 2
        for ii = 1:6
            U(NNB(ii+6, 1), ii) = U(NNB(ii+6, 1), ii) + 1;
            K(NNB(ii, 1), NNB(ii+6, 1)) = K(NNB(ii, 1), NNB(ii+6, 1)) + 1;
        end
    end
end

%% enforce Dirichlet boundary conditions if needed
if diagonal == true
    Nf = Ntot;                   % degree of freedom
    NNBc = NNB;
    dof = 1:Ntot;
else
    if problem == 0 | problem == 1 | problem == 2
        Nf = Ntot-6;
        NNBc = NNB;
        NNBc(7:12, Nx) = 0;
        dof = 1:Ntot-6;
    end
    
    % NNBc = NNB-6; % with the Dirichlet BC
    % NNBc(1:6, 1)  = 1:6;
    % NNBc(7:12, 1) = 0;

    % if problem == 0 | problem == 1 | problem == 2
    %     dof = 1:Ntot-6;     % degrees of fredom
    % elseif problem == 3
    %     dof = 7:Ntot;     % degrees of fredom
    % end
end


M = M(dof, dof); K = K(dof, dof); U = U(dof, :); 
for pp = 1:Ni
    for ee = 1:Ne
        temp1 = cell2mat(P1(pp, ee)); P1(pp, ee) = {temp1(dof, dof)};
        temp2 = cell2mat(P2(pp, ee)); P2(pp, ee) = {temp2(dof, dof)};
        temp3 = cell2mat(P3(pp, ee)); P3(pp, ee) = {temp3(dof, dof)};
        temp4 = cell2mat(Pdagger1(pp, ee)); Pdagger1(pp, ee) = {temp4(dof, dof)};
        temp5 = cell2mat(Pdagger2(pp, ee)); Pdagger2(pp, ee) = {temp5(dof, dof)};
        temp6 = cell2mat(Pdagger3(pp, ee)); Pdagger3(pp, ee) = {temp6(dof, dof)};
    end
end

%% Initial data
disp('Building the initial data..')
tic

%%% ------------------- Y0 ---------------------- %%%
Y0 = zeros(Ntot, 1);

%%% meaningless initial data %%%
% x0 = 0.5*ell;   % center
% a = 0.2*ell;    % width
% c0 = 0.1;       % magnitude
% f = @(x) c0 * exp( -1/a^2*((x-x0).^2) ) ; % bump function
% Y0 = zeros(Ntot, 1);
% for ii = 1:6
%     for kk = 1 : Nx
%         Y0(NNB(ii, kk), 1) = f(x(kk));
%     end
% end
%%% --------------------------------------------- %%%


%%% ----------- (p0, R0) and (pD, RD) ----------- %%%
p_ref = x.*[1; 0; 0];

if problem == 0 | problem == 1  
    cosa = 6/10; sina = 8/10;
    p0 = [[-cosa, -sina, 0]; [sina, -cosa, 0]; [0, 0, 1]]*(p_ref)  + [6; 0; 0];
    R0 = zeros(3, 3, Nx);
    for kk=1:Nx
        R0(:, :, kk) = [[-cosa, -sina, 0]; [sina, -cosa, 0]; [0, 0, 1]];
    end
elseif problem == 2
    p0 = p_ref;
    R0 = zeros(3, 3, Nx);
    for kk=1:Nx
        R0(:, :, kk) = eye(3);
    end
    
    RR = @(theta) [[cos(theta), -sin(theta), 0];...
                [sin(theta), cos(theta), 0];...
                [0, 0, 1]];
    RD = zeros(3, 3, Nt);
    pD = zeros(3, Nt);
    for nn = 1:Nt
        RD(:, :, nn) = RR(func_theta(t(nn)));
    end      
end

%%% test of the function func_theta
% theta_val = zeros(1, Nt);
% for nn = 1:Nt
%     theta_val(1, nn) = func_theta(t(nn));
% end
% figure();
% plot(t, theta_val);

%%% test of the function rotm2quat
%     q0_true = zeros(4, Nx); % test
%     for kk=1:Nx
%         q0_true(:, kk) = rotm2quat(R0(:, :, kk));
%     end

%%% test if we can recover p0, R0 from z0 via transfo :
% test_transfoIni(p0, R0, Y0, NNB, x, kap)

toc

%% Initialization of the state Y
Y_exp = zeros(Nf, Nt);     % zero matrix
H = zeros(4, Nt, Nx);      % quaternions
R = zeros(3, 3, Nt, Nx);   % rotation matrices
p = zeros(3, Nt, Nx);      % positions
Y_diag = zeros(Ntot, Nt);
Y_phys = zeros(Ntot, Nt);  
Bxt = zeros(3, Nt, Nx);

% initial conditions at t=0 for Y
Y_phys(:, 1) = Y0(:, 1);

if diagonal == true
    for ss = 1:Nx 
        Y_diag(NNB(:, ss), 1)= d_L * Y_phys(NNB(:, ss), 1);
    end
    Y_exp(:, 1) = Y_diag(dof, 1);  % set the initial data at time t1
else
    Y_exp(:, 1) = Y_phys(dof, 1); 
end

% initial conditions at t=0 for the position variables
for ss = 1:Nx                                 % initial data at t=0
    H(:, 1, ss) = rotm2quat(R0(:, :, ss))';   % quat angle at t=0
    R(:, :, 1, ss) = R0(:, :, ss);            % angle at t=0
    p(:, 1, ss) = p0(:, ss);                  % position at t=0
    Bxt(:, 1, ss) = R(:, :, 1, ss)*( Y_phys(NNB(1:3, ss), 1) );
end

% position variables for alternative solving centerline method:
R_ini = zeros(3, 3, Nt);   % rotation matrices
p_ini = zeros(3, Nt);      % positions
Bxt_ini = zeros(3, Nt);
R_ini(:, :, 1) = R0(:, :, 1);            % angle at t=0
p_ini(:, 1) = p0(:, 1);                  % position at t=0
Bxt_ini(:, 1) = R0(:, :, 1)*( Y_phys(NNB(1:3, 1), 1) );

%% Solving the ODE via implicit midpoint rule

% for old version of linearized solve, see simu_GEB_DIAG_V4
if linearized == true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        linearized IGEB           %
    % Approx equation: M y_t + K y = 0 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Solving the linearized system..')
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                  IGEB                    %
    % Approx equation: M y_t + K y + Q(y)y = 0 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Solving the semilinear system..')
end

%%% M.A. code %%%
tolzero=1e-12;  %-12
reltolX=1e-6;   %-6
tolF=1e-13;     %-15
%%%%%%%%%%%%%%%%%

tic

if diagonal == true
    Pi1 = 1/sqrt(2)*[zeros(3), eye(3), zeros(3), eye(3)];
else
    Pi1 = [zeros(3), eye(3), zeros(3), zeros(3)];
end
Pi2 = zeros(12, Nf);
Pi2(:, NNBc(:, 1)) = eye(12);   
Proj_v2 = Pi1*Pi2;

for kk=1:Nt-1  % loop over time

    %%% computing Y at time k+1
    Yk = Y_exp(:, kk);
    Hk = H(:, kk, 1);

    if linearized == true
        Q_yk = zeros(Nf); 
        Qdagger_yk = zeros(Nf);
    else
        [Q_yk, Qdagger_yk] = Q(Yk, P1, P2, P3, Pdagger1, Pdagger2,...
            Pdagger3, Ni, Ne, Nf, NNBc, diagonal);
    end

    zm = Yk; % initialize zm
    if approx_rot == 0 | approx_rot == 1
        xm = Hk;
        zxm = [zm; xm];
    else
        zxm = zm;
    end
%     zm = rand(Nf, 1);
%     xm = rand(4, 1);
%     zxm = [zm; xm];

    Rk_tr = transpose(func_quat2rotm(Hk));

    fk = f_ext(t(kk), problem);
    fk1 = f_ext(t(kk+1), problem);
    fk1a = fk1(1:3, 1); 
    fk1b = fk1(4:6, 1);
    JW12m_a = ht/2*U*[[Wdagger(0, fk1a), Wdagger(1, fk1a),...
        Wdagger(2, fk1a), Wdagger(3, fk1a)];...
        [Wdagger(0, fk1b), Wdagger(1, fk1b),...
        Wdagger(2, fk1b), Wdagger(3, fk1b)]];

    while 1 
        %%% Newton method: the scheme reads %%%
        % zm1 = zm - (Jac Wk(zm))^{-1} Wk(zm) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ( we call here Wm = Wk(zm) )

        if linearized == true
            Q_zm = zeros(Nf);          
            Qdagger_zm = zeros(Nf);
        else
            [Q_zm, Qdagger_zm] = Q(zm, P1, P2, P3, Pdagger1, Pdagger2,...
                Pdagger3, Ni, Ne, Nf, NNBc, diagonal);
        end

        if approx_rot == 0 | approx_rot == 1
            Rm_tr = transpose(func_quat2rotm(xm));            
            if problem == 0 | problem == 1
                    Fk1 = blkdiag(Rm_tr, Rm_tr)*fk1;
                    Fk =  blkdiag(Rk_tr, Rk_tr)*fk;           
            elseif problem == 2
                Fk1 = fk1;
                Fk =  fk;
            end            
            if approx_rot == 0
                W1m = (M + ht/2*K)*zm - (M - ht/2*K)*Yk +...
                    ht/4*(Q_yk*Yk + (Q_yk + Qdagger_yk)*zm + Q_zm*zm) +...
                    ht/2*U*( Fk1 + Fk );
            elseif approx_rot == 1
                W1m = (M + ht/2*K)*zm - (M - ht/2*K)*Yk +...
                    ht/4*(Q_yk*Yk + (Q_yk + Qdagger_yk)*zm + Q_zm*zm) +...
                    ht*U*Fk;
            end

            W2m = ( eye(4) - ht/4*func_U( Proj_v2*(zm+Yk) ))*xm -...
                ( eye(4) + ht/4*func_U( Proj_v2*(zm+Yk) ) )*Hk;            

            Wm = [W1m; W2m];

            JW11m = M + ht/2*K + ht/4*(Q_yk + Qdagger_yk) +...
                ht/4*(Q_zm + Qdagger_zm);
            if problem == 0 | problem == 1
                JW12m = JW12m_a*blkdiag(xm, xm, xm, xm);
            elseif problem == 2 | approx_rot == 1
                JW12m = zeros(Nf, 4);
            end
            JW21m = -ht/4*func_Udagger( xm + Hk )*Proj_v2;
            JW22m = eye(4) - ht/4*func_U( Proj_v2*(zm+Yk) );            
            JWm = [[JW11m, JW12m]; [JW21m, JW22m]];
        else
            if problem == 0 | problem == 1
                Fk =  blkdiag(Rk_tr, Rk_tr)*fk;           
            elseif problem == 2
                Fk =  fk;
            end
            Wm = (M + ht/2*K)*zm - (M - ht/2*K)*Yk +...
                    ht/4*(Q_yk*Yk + (Q_yk + Qdagger_yk)*zm + Q_zm*zm) +...
                    ht*U*Fk;
            JWm = M + ht/2*K + ht/4*(Q_yk + Qdagger_yk) +...
                ht/4*(Q_zm + Qdagger_zm);
        end
        
        zxm1 = zxm - JWm\Wm; 

        %%% M.A. code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rel_err = (zxm1 - zxm)./zxm;                                     %
        nan_or_inf = find( isnan(rel_err) + isinf(rel_err) ...        %
            + (abs(zxm)<=tolzero) );                                   %
        rel_err(nan_or_inf) = 0;                                      %
        if (norm(rel_err,inf) <= reltolX) && (norm(Wm,inf) <= tolF)%
            break                                                     %
        end                                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

        zxm = zxm1;
        zm = zxm1(1:Nf, 1);
        if approx_rot == 0 | approx_rot == 1
            xm = zxm1(Nf+1:Nf+4, 1);
        end
    end
    Y_exp(:, kk+1) = zm; % we just computed Y at time k+1
    
    if diagonal == true
        Y_diag(dof, kk+1) = Y_exp(:, kk+1);
        Y_phys(NNB(:, 1), kk+1) = d_Linv*Y_diag(NNB(:, 1), kk+1); % at x=0
    else
        Y_phys(dof, kk+1) = Y_exp(:, kk+1);
    end
    
    if approx_rot == 0 | approx_rot == 1
        H(:, kk+1, 1) = xm;
        R_ini(:, :, kk+1) = func_quat2rotm(xm);

        % get the position information at x = 0
        idx = 1:kk+1;
        v1_k1 = Y_phys(NNB(1:3, 1), kk+1);
        Bxt_ini(:, kk+1) = R_ini(:, :, kk+1)*v1_k1;
        p_ini(1, kk+1) = p_ini(1, 1) + trapz(t(idx), Bxt_ini(1, idx), 2);
        p_ini(2, kk+1) = p_ini(2, 1) + trapz(t(idx), Bxt_ini(2, idx), 2);
        p_ini(3, kk+1) = p_ini(3, 1) + trapz(t(idx), Bxt_ini(3, idx), 2);
    end

    %%% now we recover the position of the beam at this time k+1 %%%
    for ss = 1:Nx   % for all x    
        % recover physical variable
        if diagonal == true
            Y_phys(NNB(:, ss), kk+1) = d_Linv*Y_diag(NNB(:, ss), kk+1);
        end
        
        %%% orientation of the cross sections %%%
        v2_k = Y_phys(NNB(4:6, ss), kk);
        v2_k1 = Y_phys(NNB(4:6, ss), kk+1);
        Axt = func_U( ( v2_k + v2_k1 )/2);  % uses angular velocities  
        H(:, kk+1, ss) = (eye(4)-ht/2*Axt)\((eye(4) + ht/2*Axt)*H(:,kk, ss)); % GOOD
        R(:, :, kk+1, ss) = func_quat2rotm(H(:, kk+1, ss));

        %%% position of the reference line %%%
        idx = 1:kk+1;
        v1_k1 = Y_phys(NNB(1:3, ss), kk+1);
        Bxt(:, kk+1, ss) = R(:, :, kk+1, ss)*v1_k1;
        p(1, kk+1, ss) = p(1, 1, ss) + trapz(t(idx), Bxt(1, idx, ss), 2);
        p(2, kk+1, ss) = p(2, 1, ss) + trapz(t(idx), Bxt(2, idx, ss), 2);
        p(3, kk+1, ss) = p(3, 1, ss) + trapz(t(idx), Bxt(3, idx, ss), 2);
    end
    %%% ------------------------------------------------------ %%% 
end
toc

%% Full state including the Dirichlet nodes
%ymax = max(max(Y)); ymin = min(min(Y));

%% True approximation of y ?
% % % disp('Building the approximation of y..')
% % % tic
% % % 
% % % %xNew = linspace(0, ell, 2000); % different spatial discretization
% % % Nx_true = Nx*50;
% % % xNew = linspace(0, ell, Nx_true); % same spatial discretization
% % % y_true = zeros(Nx_true*Ni, Nt);
% % % NNB_true = reshape(1:Nx_true*Ni, Ni, Nx_true);
% % % for ii = 1:Ni
% % %     for nn = 1:Nt
% % %         for kk = 1:numel(xNew)
% % %             y_true(NNB_true(ii, kk), nn) = func_yi(xNew(kk), ii, nn, Y_phys, Nx, NNB,Ne,he, x);
% % %         end
% % %     end
% % % end
% % % p_ref_true = xNew.*[1; 0; 0];
% % % 
% % % if problem == 0 | problem == 1  
% % %     cosa = 6/10; sina = 8/10;
% % %     p0_true = [[-cosa, -sina, 0]; [sina, -cosa, 0]; [0, 0, 1]]*(p_ref_true)  + [6; 0; 0];
% % %     R0_true = zeros(3, 3, Nx_true);
% % %     for kk=1:Nx_true
% % %         R0_true(:, :, kk) = [[-cosa, -sina, 0]; [sina, -cosa, 0]; [0, 0, 1]];
% % %     end
% % % % elseif problem == 2
% % % %     theta = 1.5;
% % % %     RR = [[cos(theta), -sin(theta), 0];...
% % % %         [sin(theta), cos(theta), 0];...
% % % %         [0, 0, 1]];
% % % % 
% % % %     p0 = p_ref;
% % % %     pD = RR*p_ref;
% % % %     
% % % %     R0 = zeros(3, 3, Nx);
% % % %     RD = zeros(3, 3, Nx);
% % % %     for kk=1:Nx
% % % %         R0(:, :, kk) = eye(3);
% % % %         RD(:, :, kk) = RR;
% % % %     end
% % % end
% % % 
% % % % truc = 0;
% % % % for ii=1:12
% % % %     truc = truc + sum(sum(Y(NNB(ii, :), :)-y(:, :, ii)));
% % % % end
% % % % disp(['truc = ', num2str(truc)])
% % % 
% % % toc

%% Plot the solution
if plot_y
    disp('Plots of the solution y physical..')
    f = figure();
    set(gcf,'Position',[100 100 1200 600])
    ii2subplot = [1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12];
    ii2title = ["Linear velocity", "Angular velocity", "Forces", "Moments"];
    ii2label = ["$V_1$","$V_2$","$V_3$","$W_1$","$W_2$","$W_3$",...
        "$\Phi_1$","$\Phi_2$","$\Phi_3$","$\Psi_1$","$\Psi_2$","$\Psi_3$"];
    for ii = 1:12
        subplot(3, 4, ii2subplot(ii));
        s = surf(x, t, Y_phys(NNB(ii, :), :)');
        s.EdgeColor = 'none'; 
        ylabel([ii2label(ii),'\ \ $t$'],'Interpreter','latex'); 
        xlabel('$x$','Interpreter','latex'); 
        %zlabel(ii2label(ii),'Interpreter','latex');
        colorbar
        view(2)       
        if ii2subplot(ii) <= 4
            title(ii2title(ii2subplot(ii)),'Interpreter','latex', 'fontsize', 12);
        end
    end
    orient(f,'landscape')
%     exportgraphics(f,'flying_spaghetti_Y_view2.pdf','ContentType','vector')

    if diagonal == true
        disp('Plots of the solution y diagonal..')
        f = figure();
        set(gcf,'Position',[100 100 1200 600])
        ii2subplot = [1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12];
        ii2title = ["r-1", "r-2", "r+1", "r+2"];
        ii2label = ["$V_1$","$V_2$","$V_3$","$W_1$","$W_2$","$W_3$",...
            "$\Phi_1$","$\Phi_2$","$\Phi_3$","$\Psi_1$","$\Psi_2$","$\Psi_3$"];
        for ii = 1:12
            subplot(3, 4, ii2subplot(ii));
            s = surf(x, t, Y_diag(NNB(ii, :), :)');
            s.EdgeColor = 'none'; 
            ylabel([ii2label(ii),'\ \ $t$'],'Interpreter','latex'); 
            xlabel('$x$','Interpreter','latex'); 
            %zlabel(ii2label(ii),'Interpreter','latex');
            colorbar
            view(2)
            if ii2subplot(ii) <= 4
                title(ii2title(ii2subplot(ii)),'Interpreter','latex', 'fontsize', 12);
            end
        end
        orient(f,'landscape')
    %     exportgraphics(f,'flying_spaghetti_Y_view2.pdf','ContentType','vector')
    end
end

%% H1 norm of y
if plot_norm == true
    Y_H0 = zeros(1, Nt);
    Y_H1 = zeros(1, Nt);
    for kk = 1:Nt
        Y_H0(1, kk) = (Y_phys(:, kk)')*LL0*Y_phys(:, kk);
        Y_H1(1, kk) = (Y_phys(:, kk)')*LL0*Y_phys(:, kk) + (Y_phys(:, kk)')*LL1*Y_phys(:, kk);
    end
    figure()
    plot(t, Y_H0);
    hold on;
    plot(t, Y_H1);
    legend('L2 norm of y', 'H1 norm of y');
end

%% Y at boundary
if problem == 0 | problem == 1
    Z1 = zeros(3, Nt);
    Z2 = zeros(3, Nt);
    for kk = 1:Nt
        Z1(:, kk) = R(:, :, kk, 1)*Y_phys(NNB(7:9, 1), kk);
        Z2(:, kk) = R(:, :, kk, 1)*Y_phys(NNB(10:12, 1), kk);
    end
    figure()
    plot(t, Z1(1, :));
    hold on;
    plot(t, Z2(2, :));
    plot(t, Z2(3, :));
    legend('7, x=0', '8, x=0', '12, x=0');
    title('Value of Y at x=0')
elseif problem == 2
    figure()
    plot(t, Y_phys(NNB(6, 1), :));
    legend('6, x=0');
    title('Value of Y at x=0')
end

%% test other way recover centerline
if true
    centerline_mode = "XSolve";      % or "XSolve"
%     [p2, R2] = recover_position(p_ini(:, :), R_ini(:, :, :), Y_phys, NNB, centerline_x_scheme, centerline_mode, x, t, flexMat, kap);
    [p2, R2] = recover_position(p(:, :, 1), R(:, :, :, 1), Y_phys, NNB, centerline_x_scheme, centerline_mode, x, t, flexMat, kap);
%     [p2, R2] = recover_position(pD, RD, Y_phys, NNB, centerline_x_scheme, centerline_mode, x, t, flexMat, kap);
end

%% permute space and time
p = permute(p, [1, 3, 2]); % for plotting we need to change the order of the space and time indexes

%% position of the beam through time
% centerline_mode = "TSolve";      % or "XSolve"
% [p, R] = recover_position(p0, R0, Y_phys, NNB, centerline_t_scheme, centerline_mode, x, t, flexMat, kap);
% % % [p, R] = recover_position(p0_true, R0_true, y_true, NNB_true, centerline_t_scheme, centerline_mode, xNew, t, flexMat, kap);


%%% ----------- plot arclength ----------- %%%
file_name_arclength = ['fig/ARCLEN_p_', linNonlin, '_', centerline_mode, '.pdf'];
plot_arclen(p, x, t, centerline_mode, file_name_arclength);
file_name_arclength = ['fig/ARCLEN_p2_', linNonlin, '_', centerline_mode, '.pdf'];
plot_arclen(p2, x, t, centerline_mode, file_name_arclength);


% arcLength_time = zeros(Nt, 1);
% for nn = 1:Nt
%     arcLength_time(nn, 1) = trapz(x, sqrt(sum((gradient(p(:, :, nn), hx)).^2))); % arclength
% end
% f_arc = figure();
% plot(t, arcLength_time, 'lineWidth', 2);
% xlabel('$t$','Interpreter','latex');
% grid on
% title(title_arclength, 'Interpreter', 'latex');
% %print(f_arc, file_name_arclength,'-dpdf')
% %exportgraphics(f_arc,file_name_arclength,'ContentType','vector');
%%% --------------------------------------- %%%

%%% ----------- plot centerline ----------- %%%
if centerline_mode == 'XSolve'
    title_centerline = 'Centerline: by space integration using $z$';
elseif centerline_mode == 'TSolve'
    title_centerline = 'Centerline: by time integration using $v$';
end
file_name_centerline = ['CENTERL_', linNonlin, '_', centerline_mode, '.pdf'];


%pf = zeros(3, Nx); % the undeformed configuration fulfilling the clamped BC
% temp = R0(:, 1, 1);
% for kk = 1:Nx
%     pf(:, kk) = p0(:, 1) + x(kk)*temp;
% end
if ht == 0.01
    fact = 10;
elseif ht == 0.005
	fact = 20;
else
    fact = 1;
end
% fact = 1
%[-37.5, 30];
if plot_centerline
    if problem == 0
        centerline_t = 1:(5*fact):Nt;  % time at which we plot the centerline
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        titleCenterline = 'Flying spaghetti problem 2D';
        f_c = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        centerline_t = 1:(5*fact):Nt;  % time at which we plot the centerline
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        titleCenterline = 'Flying spaghetti problem 2D p2';
        f_c = plotCenterline(p2, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        
    elseif problem == 1
        titleCenterline = 'Flying spaghetti problem 3D';
        
%         centerline_t = fix([0, 2, 3, 3.8, 4.4, 5, 5.5, 5.8, 6.1, 6.5]*10*fact+1);
%         viewCent = [0, 90];
%         %viewCent = [90, 0];
%         locLegend = 'northeastoutside';
%         f_c = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx);
% 
%         centerline_t = fix([0, 2.5, 3.5, 3.8, 4.5]*100+1);
%         viewCent = [90, 0];
%         locLegend = 'southwest';
%         f_c = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx);
%         camroll(90);
        
        centerline_t = 1:(4*fact):Nt;
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
    
        titleCenterline = 'Flying spaghetti problem 3D p2';
        centerline_t = 1:(4*fact):Nt;
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c = plotCenterline(p2, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
    
    elseif problem == 2
        titleCenterline = 'Book p fin';
        centerline_t = (5*10*fact+1):(5*fact):(9*10*fact+1);
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c1 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        titleCenterline = 'Book p debut';
        centerline_t = 1:(5*fact):(5*10*fact+1);
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c2 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        
        titleCenterline = 'Book p2 fin';
        centerline_t = (5*10*fact+1):(5*fact):(9*10*fact+1);
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c3 = plotCenterline(p2, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        titleCenterline = 'Book p2 debut';
        centerline_t = 1:(5*fact):(5*10*fact+1);
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c4 = plotCenterline(p2, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
%         centerline_t = (5.5*10*fact+1):(5*fact):(9*10*fact+1);
%         viewCent = [0, 90];
%         locLegend = 'northeastoutside';
%         f_c = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
    end

%     f_c = figure();
%     h = plot3(p0(1, :), p0(2, :), p0(3, :), 'g', 'lineWidth', 2); 
%     hold on;
%     %h.Color(4) = 0.5;
%     grid on;
%     sec = centerline_t(2);
%     h = plot3(p(1, :, sec), p(2, :, sec), p(3, :, sec), 'k', 'lineWidth', 1);
%     for nn = centerline_t(3:end-1)
%         h = plot3(p(1, :, nn), p(2, :, nn), p(3, :, nn), 'k', 'lineWidth', 1);
%         set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     end
%     plot3(p(1, :, end), p(2, :, end), p(3, :, end), 'r', 'lineWidth', 2);
%     %plot3(pf(1, :), pf(2, :), pf(3, :), '--g', 'lineWidth', 3);
% 
%     pp = permute(p, [1, 3, 2]);
%     plot3(pp(1, :, 1), pp(2, :, 1), pp(3, :, 1), '--b');
%     plot3(pp(1, :, Nx), pp(2, :, Nx), pp(3, :, Nx), '--r');
% 
%     axis equal;
%     view(gca, orient_cent(ii, 1), orient_cent(ii, 2));
%     legend('$t=0$', '$t$', ['$t = ', num2str(T), '$'],...
%         '$x=0$', '$x=\ell$','Interpreter','latex',...
%         'Location', locLegend(ii), 'fontsize', 12);
%     xlabel('X', 'Interpreter', 'latex');
%     ylabel('Y', 'Interpreter', 'latex');
%     zlabel('Z', 'Interpreter', 'latex');
%     title('Flying spaghetti problem', 'Interpreter', 'latex', 'fontsize', 12);
%     %exportgraphics(f_c,'flying_spaghetti.pdf','ContentType','vector')
end
%%% -------------------------------------- %%%
disp('End.')




function [res, resDagger] = Q(Y, P1, P2, P3, Pdagger1, Pdagger2, Pdagger3, Ni, Ne, Nf, NNB, diagonal)
% Q computes both Q and Q_dagger
% both are of size Nf x Nf
    res=zeros(Nf); resDagger=zeros(Nf);

    if diagonal == true
        for pp=1:Ni
            for ee=1:Ne
                res = res + Y(NNB(pp, 2*ee-1))*cell2mat(P1(pp, ee)) ...
                          + Y(NNB(pp, 2*ee))*cell2mat(P2(pp, ee)) ...
                          + Y(NNB(pp, 2*ee+1))*cell2mat(P3(pp, ee));
                resDagger = resDagger + Y(NNB(pp, 2*ee-1))*cell2mat(Pdagger1(pp, ee)) ...
                                      + Y(NNB(pp, 2*ee))*cell2mat(Pdagger2(pp, ee)) ...
                                      + Y(NNB(pp, 2*ee+1))*cell2mat(Pdagger3(pp, ee));
            end
        end
    else
        for pp=1:6
            for ee=1:Ne
                res = res + Y(NNB(pp, 2*ee-1))*cell2mat(P1(pp, ee)) ...
                          + Y(NNB(pp, 2*ee))*cell2mat(P2(pp, ee)) ...
                          + Y(NNB(pp, 2*ee+1))*cell2mat(P3(pp, ee));
                resDagger = resDagger + Y(NNB(pp, 2*ee-1))*cell2mat(Pdagger1(pp, ee)) ...
                                      + Y(NNB(pp, 2*ee))*cell2mat(Pdagger2(pp, ee)) ...
                                      + Y(NNB(pp, 2*ee+1))*cell2mat(Pdagger3(pp, ee));
            end
        end
        for pp=7:Ni
            for ee = 1:Ne-1
                res = res + Y(NNB(pp, 2*ee-1))*cell2mat(P1(pp, ee)) ...
                          + Y(NNB(pp, 2*ee))*cell2mat(P2(pp, ee)) ...
                          + Y(NNB(pp, 2*ee+1))*cell2mat(P3(pp, ee));
                resDagger = resDagger + Y(NNB(pp, 2*ee-1))*cell2mat(Pdagger1(pp, ee)) ...
                                      + Y(NNB(pp, 2*ee))*cell2mat(Pdagger2(pp, ee)) ...
                                      + Y(NNB(pp, 2*ee+1))*cell2mat(Pdagger3(pp, ee));
            end
            ee=Ne;
            res = res + Y(NNB(pp, 2*ee-1))*cell2mat(P1(pp, ee)) ...
                      + Y(NNB(pp, 2*ee))*cell2mat(P2(pp, ee));
            resDagger = resDagger + Y(NNB(pp, 2*ee-1))*cell2mat(Pdagger1(pp, ee)) ...
                                  + Y(NNB(pp, 2*ee))*cell2mat(Pdagger2(pp, ee));
        end
    end
end