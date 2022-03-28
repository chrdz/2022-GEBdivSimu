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
clear
close all
clc
fs = 14; % font size
size_line = 2; % line size for plots;

%% To be chosen
% Do you want the linearized system?
linearized = false;  % true: we take the IGEB system without g(y)
                     % false: we take the real (semilinear) IGEB system
            
% What do you want to plot?
plot_y = true;          % true: we plot the solution y
                        % false: we don't plot it
plot_centerline = true; % true: we plot the centerline of the beam
                        % false: we don't plot it                        
plot_norm = false;      % true: we plot the norm of y(., t) through time
                        % we don't plot it
plot_boundary = false;   % we plot the value of y at x=0 through time
                        % we don't plot it

% Which scheme for recovering displacement variables?                        
centerline_t_scheme = 0;     % 0: mid point rule
                             % 1: explicit Euler (does not work properly)
                             % 2: ode45 (RK)
                             % 3: implicit Euler
                             % 4: zupan paper scheme (does not work properly)           
centerline_x_scheme = 0;     % 0: mid point rule
                             % 1: explicit Euler (does not work properly)

% Which problem do you want to solve
problem = 3;    % 0: flying spaghetti problem 2D // i.e. Problem II - 2D
                % 1: flying spaghetti problem 3D // i.e. Problem II - 3D
                % 2: rotating arm problem // i.e. Problem I
                % 3: feedback control problem // i.e. Problem III
                
BC = 1;   % 0: we take transparent boundary conditions (BC)
          % 1: we take something close to transparent BC
          % 2: we take somthing far from transparent BC   (does not work)
          
zeroOrderBC = false; % true: we choose the initial velocities so that the 
                    % the zero-order compatibility conditions of the IGEB
                    % model are fulfilled
                    % false: the initial velocities are set to zero

% How do you write the Newton scheme?
approx_rot = 0;    % 0 most precise: mid point everywhere
                   % 1 no mid point for Rfext in W1m
                   % 2 no mid point for Rfext in W1m and in calU(v2)
if problem == 3
    approx_rot = 2;
end

% Do you want to work with the diagonal system?      
diagonal = false;   % true: we diagonalize the system before solving
                   % false: we work directly with the physical system

if problem == 0 || problem == 1 || problem == 2                   
    ell = 10;    % length of the space interval
    Ne = 20;     % number of elements
    T = 15;      % end of the time interval
    ht = 0.1;   % time step
elseif problem == 3
    ell = 1;
    Ne = 15;
    Tfactor = 5;
    T = 1;         % end of the time interval
    ht = 0.002;    % time step
end
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
if problem == 0 || problem == 1 || problem == 3
    EI = 500; GJ = 500;
elseif problem == 2
    EI = 1000; GJ = 1000;
end
rhoA = 1; rhoJ = diag([20, 10, 10]);
massMat = blkdiag(rhoA*eye(3), rhoJ);            % the MASS matrix
flexMat = inv(diag([EA, GAs, GAs, GJ, EI, EI])); % the FLEXIBILITY matrix

if problem == 3 % build feedback matrix K
    K_fb = eye(6)*flexMat^(-1/2)*massMat^(1/2); % transparent BC
    diag_kappa = diag(K_fb);
    mu1 = sqrt(max(diag_kappa(1:3)*min(diag_kappa(1:3))));
    mu2 = sqrt(max(diag_kappa(4:6)*min(diag_kappa(4:6))));
    new_kappa_diag = [mu1, mu1, mu1, mu2, mu2, mu2];
    if BC == 1          % close to transparent BC
        K_fb = diag(new_kappa_diag);
    elseif BC == 2      % far from transparent BC
        K_fb = K_fb +  diag(rand(1, 6))*diag([0.001, 0.001, 0.001, 0.1, 0.1, 0.1])*5;
    end
end

%% Coeffficients of the PDE system
J = blkdiag(massMat, flexMat);
B = zeros(12, 12); B(1:6, 7:12) = - boldE; B(7:12, 1:6) = transpose(boldE);
A = zeros(12, 12); A(1:6, 7:12) = -eye(6); A(7:12, 1:6) = -eye(6);
Ni = 12;    % number of PDEs in our system

PiMinus = [eye(6), zeros(6)];
PiPlus = [zeros(6), eye(6)];
if diagonal == false
    if problem == 2 % rotating arm
        Wa = (PiMinus') * PiPlus;
        Wb = PiPlus';
    elseif problem == 0 || problem == 1 % flying spaghetti
        Wa = (PiPlus') * PiMinus;
        Wb = -PiMinus';
    elseif problem == 3
        Wa = (PiMinus') * K_fb * PiMinus - (PiPlus') * PiMinus;
    end
end

%% Coeffficients of the PDE system in DIAGONAL form
if diagonal == true
    d_L = 1/sqrt(2)*[[eye(6), eye(6)]; [eye(6), -eye(6)]];
    d_Linv = d_L;
    d_J = d_L*J*d_Linv;
    d_D = [[-eye(6), zeros(6)]; [zeros(6), eye(6)]];
    d_G = @(r, M, C) d_L*G(d_Linv*r, M, C)*d_Linv;
    d_Gdagger = @(r, M, C) d_L*Gdagger(d_Linv*r, M, C)*d_Linv;
    d_B = d_L*B*d_Linv;
    if problem == 0 || problem == 1
        d_Wa = d_D*[[zeros(6), eye(6)]; [zeros(6), eye(6)]];
        d_Wb = -d_D*[[eye(6), zeros(6)]; [eye(6), zeros(6)]];
        d_Wc = -sqrt(2)*[zeros(6); eye(6)];
    elseif problem ==  2
        d_Wa = d_D*[[zeros(6), eye(6)]; [zeros(6), eye(6)]];
        d_Wb = -d_D*[[eye(6), zeros(6)]; [-eye(6), zeros(6)]];
        d_Wc = -sqrt(2)*[zeros(6); eye(6)];
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
Z = sparse(Ntot, 6);

% % % if diagonal == true
% % %     for ii = 1:Ni
% % %         for jj = 1:Ni
% % %             K(NNB(ii, Nx), NNB(jj, Nx)) = K(NNB(ii, Nx), NNB(jj, Nx)) + d_Wa(ii, jj);
% % %             K(NNB(ii, 1), NNB(jj, 1)) = K(NNB(ii, 1), NNB(jj, 1)) + d_Wb(ii, jj);
% % %         end
% % %         for jj = 1:6
% % %             Z(NNB(ii, 1), jj) = Z(NNB(ii, 1), jj) + d_Wc(ii, jj);
% % %         end
% % %     end
% % % else
% % %     if problem == 0 || problem == 1
% % %         for ii = 1:6
% % %             Z(NNB(ii, 1), ii) = Z(NNB(ii, 1), ii) - 1;
% % %             K(NNB(ii+6, 1), NNB(ii, 1)) = K(NNB(ii+6, 1), NNB(ii, 1)) + 1;
% % %         end
% % %     elseif problem == 2
% % %         for ii = 1:6
% % %             Z(NNB(ii+6, 1), ii) = Z(NNB(ii+6, 1), ii) + 1;
% % %             K(NNB(ii, 1), NNB(ii+6, 1)) = K(NNB(ii, 1), NNB(ii+6, 1)) + 1;
% % %         end
% % %     end
% % % end

if diagonal == true
    for ii = 1:Ni
        for jj = 1:Ni
            K(NNB(ii, Nx), NNB(jj, Nx)) = K(NNB(ii, Nx), NNB(jj, Nx)) + d_Wa(ii, jj);
            K(NNB(ii, 1), NNB(jj, 1)) = K(NNB(ii, 1), NNB(jj, 1)) + d_Wb(ii, jj);
        end
        for jj = 1:6
            Z(NNB(ii, 1), jj) = Z(NNB(ii, 1), jj) + d_Wc(ii, jj);
        end
    end
else
    if problem == 0 || problem == 1 || problem == 2
        for ii = 1:Ni
            for jj = 1:Ni
                K(NNB(ii, 1), NNB(jj, 1)) = K(NNB(ii, 1), NNB(jj, 1)) + Wa(ii, jj);
            end
            for jj = 1:6
                Z(NNB(ii, 1), jj) = Z(NNB(ii, 1), jj) + Wb(ii, jj);
            end
        end
    elseif problem == 3
        for ii = 1:Ni
            for jj = 1:Ni
                K(NNB(ii, Nx), NNB(jj, Nx)) = K(NNB(ii, Nx), NNB(jj, Nx)) + Wa(ii, jj);
            end
        end
    end
end
    

%% enforce Dirichlet boundary conditions if needed
if diagonal == true
    Nf = Ntot;                   % degree of freedom
    NNBc = NNB;
    dof = 1:Ntot;
else
    if problem == 0 || problem == 1 || problem == 2
        % free beam at x=ell
        Nf = Ntot-6;
        NNBc = NNB;
        NNBc(7:12, Nx) = 0;
        dof = 1:Ntot-6;
    elseif problem == 3
        % clamped beam at x=0
        Nf = Ntot-6;
        NNBc = NNB-6; 
        NNBc = subplus(NNBc);
        dof = 7:Ntot;
    end
    % free beam at x=0
    % NNBc = NNB-6; % with the Dirichlet BC
    % NNBc(1:6, 1)  = 1:6;
    % NNBc(7:12, 1) = 0;
end


M = M(dof, dof); K = K(dof, dof); Z = Z(dof, :); 
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
if problem == 0 || problem == 1 || problem == 2
    Y0 = zeros(Ntot, 1);
elseif problem == 3
    Gamma0 = [0; 0; 0];    % there is no initial shear
    W0hat = [0, -1/sqrt(2), 0; 1/sqrt(2), 0, 1/sqrt(2); 0, -1/sqrt(2), 0];
    Upsilon0 = func_vec(W0hat);
    z0 = [Gamma0; Upsilon0];  % strains
    z0 = flexMat\z0;          % corresponding stresses
    v0_ell = - K_fb\z0;      % velocities at ell fulfilling 
                              % compatibility conditions at x = ell
    %%% imposes zero-order comp. cond + null acceleration at x=0:
    y0Mat = zeros(12, Nx);
    if zeroOrderBC == true
        for ii = 1:6
            x_constr = [0, 0.05, ell-0.05, ell];
            y_constr = [0, 0, v0_ell(ii), v0_ell(ii)];
            v0_interp = pchip(x_constr, y_constr, x);
            for kk = 1:Nx
                y0Mat(ii, :) = v0_interp;
            end
        end
    end
    for kk = 1:Nx
        y0Mat(7:12, kk) = z0;
    end
    %%% transform y0Mat to a vector:
    Y0 = zeros(Ntot, 1);
    for ii = 1:Ni
        for kk = 1 : Nx
            Y0(NNB(ii, kk), 1) = y0Mat(ii, kk);
        end
    end
end
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

if problem == 0 || problem == 1
    % in problems 0 and 1 the beam is free at both ends
    % hence there is no boundary data
    
    % simo version:
    cosa = 6/10; sina = 8/10;
    p0 = [[-cosa, -sina, 0]; [sina, -cosa, 0]; [0, 0, 1]]*(p_ref)  + [6; 0; 0];
    R0 = zeros(3, 3, Nx);
    for kk=1:Nx
        R0(:, :, kk) = [[-cosa, -sina, 0]; [sina, -cosa, 0]; [0, 0, 1]];
    end
    
%     % hesse version TEST:
%     cosa = 6/10; sina = 8/10;
%     p0 = [[-cosa, 0, -sina]; [0, 1, 0]; [sina, 0, -cosa]]*(p_ref);
%     R0 = zeros(3, 3, Nx);
%     for kk=1:Nx
%         R0(:, :, kk) = [[-cosa, 0, -sina]; [0, 1, 0]; [sina, 0, -cosa]];   
%     end
    
elseif problem == 2
    % in problem 2, the beam is clamped and the angle changes with time
    % hence we specify both the initial and boundary data
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
    
elseif problem == 3
    syms eta;     
    p0_syms = 1/sqrt(2)*[eta; (1-cos(eta)); sin(eta)];
    R0_syms = [[1/sqrt(2), 0, -1/sqrt(2)]; 
               [sin(eta)/sqrt(2), cos(eta), sin(eta)/sqrt(2)];
               [cos(eta)/sqrt(2), -sin(eta), cos(eta)/sqrt(2)]];
    p0 = zeros(3, Nx);
    R0 = zeros(3, 3, Nx);
    for kk=1:Nx
        p0(:, kk) = subs(p0_syms, x(kk));
        R0(:, :, kk) = subs(R0_syms, x(kk));
    end
    RD = zeros(3, 3, Nt);
    pD = zeros(3, Nt);
    for nn = 1:Nt
        pD(:, nn) = p0(:, 1);
        RD(:, :, nn) = R0(:, :, 1);
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
Y_diag = zeros(Ntot, Nt);  % diagonal variable
Y_phys = zeros(Ntot, Nt);  % physical variable
Bxt = zeros(3, Nt, Nx);    % matrix were we will store some information

% initial conditions at t=0 for the variable Y
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

% position/rotation variables for alternative solving centerline method:
R_ini = zeros(3, 3, Nt);       % rotation matrices
p_ini = zeros(3, Nt);          % positions
Bxt_ini = zeros(3, Nt);
R_ini(:, :, 1) = R0(:, :, 1);  % rotation at time t=0
p_ini(:, 1) = p0(:, 1);        % position at time t=0
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

if approx_rot == 0 || approx_rot == 1
    % Definition of the matrices Pi1 and Pi2
    % --> v_2(x,t) = Pi1 * Pi2 * y(x,t)
    if diagonal == true
        Pi1 = 1/sqrt(2)*[zeros(3), eye(3), zeros(3), eye(3)];
    else
        Pi1 = [zeros(3), eye(3), zeros(3), zeros(3)];
    end
    Pi2 = zeros(12, Nf);
    Pi2(:, NNBc(:, 1)) = eye(12);
    Proj_v2 = Pi1*Pi2;
end

for kk=1:Nt-1  % loop over time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % here we want to compute Y(t_{k+1}) 
    % and at the same time R(0, t_{k+1}) in some cases
    %-------------------------------------%
    % Newton method: the scheme reads:
    % zm1 = zm - (Jac Wk(zm))^{-1} Wk(zm) %
    %-------------------------------------%
    % Here we call:  Wm  = Wk(zm)
    %     and        JWm = Jac Wk(zm) (Jacobian matrix)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Yk = Y_exp(:, kk);   % Y(t_k)
    Hk = H(:, kk, 1);    % R(0, t_k)

    if linearized == true % if the system is linearized -> no term Q(Y)Y
        Q_yk = zeros(Nf); 
        Qdagger_yk = zeros(Nf);
    else % we compute the needed variables at time t_k for latter on
        [Q_yk, Qdagger_yk] = Q(Yk, P1, P2, P3, Pdagger1, Pdagger2,...
            Pdagger3, Ni, Ne, Nf, NNBc, diagonal, problem);
    end

    % Initializing the variable zxm for the Newton scheme below
    zm = Yk; 
    if approx_rot == 0 || approx_rot == 1
        xm = Hk; zxm = [zm; xm];
    else % i.e. approx_rot == 2
        zxm = zm;
    end
%     zm = rand(Nf, 1); xm = rand(4, 1); zxm = [zm; xm]; % test

    Rk_tr = transpose(func_quat2rotm(Hk)); % R(0, t_k)^T

    % some variables involving the external forces
    fk = f_ext(t(kk), problem);    % f(t_k)
    fk1 = f_ext(t(kk+1), problem); % f(t_{k+1})
    fk1a = fk1(1:3, 1);            % f(t_k)     first 3 components
    fk1b = fk1(4:6, 1);            % f(t_{k+1}) last 3 components
      
    JW12m_a = ht/2*Z*[[Wdagger(0, fk1a), Wdagger(1, fk1a),...
        Wdagger(2, fk1a), Wdagger(3, fk1a)];...
        [Wdagger(0, fk1b), Wdagger(1, fk1b),...
        Wdagger(2, fk1b), Wdagger(3, fk1b)]];

    % the Newton loop:
    while 1 
        if linearized == true
            Q_zm = zeros(Nf);          
            Qdagger_zm = zeros(Nf);
        else
            [Q_zm, Qdagger_zm] = Q(zm, P1, P2, P3, Pdagger1, Pdagger2,...
                Pdagger3, Ni, Ne, Nf, NNBc, diagonal, problem);
        end

        % Computing Wm
        if approx_rot == 0 || approx_rot == 1 % only for problem == 0,1 or 2
            % for problem == 0 or 1, the most precise approach (i.e. close
            % to the PDE system) is with approx_rot = 0
            Rm_tr = transpose(func_quat2rotm(xm));            
            if problem == 0 || problem == 1
                    Fk1 = blkdiag(Rm_tr, Rm_tr)*fk1;
                    Fk =  blkdiag(Rk_tr, Rk_tr)*fk;
            elseif problem == 2
                Fk1 = fk1; Fk =  fk;
            end   
            if approx_rot == 0
                Fk_approx = (Fk1 + Fk)/2;
            else  % i.e. approx_rot == 1
                Fk_approx = Fk;
            end
            
            W1m = (M + ht/2*K)*zm - (M - ht/2*K)*Yk +...
                ht/4*(Q_yk*Yk + (Q_yk + Qdagger_yk)*zm + Q_zm*zm) +...
                ht*Z*Fk_approx;
            W2m = ( eye(4) - ht/4*func_U( Proj_v2*(zm+Yk) ))*xm -...
                ( eye(4) + ht/4*func_U( Proj_v2*(zm+Yk) ) )*Hk;            
            Wm = [W1m; W2m];

            % Computing Jm
            JW11m = M + ht/2*K + ht/4*(Q_yk + Qdagger_yk) +...
                ht/4*(Q_zm + Qdagger_zm);
            if problem == 0 || problem == 1
                JW12m = JW12m_a*blkdiag(xm, xm, xm, xm);
            elseif problem == 2 || approx_rot == 1
                JW12m = zeros(Nf, 4);
            end
            JW21m = -ht/4*func_Udagger( xm + Hk )*Proj_v2;
            JW22m = eye(4) - ht/4*func_U( Proj_v2*(zm+Yk) );            
            JWm = [[JW11m, JW12m]; [JW21m, JW22m]];
            
        else  % i.e. approx_rot == 2
              % we don't look for R(0, t_k) at the same time
              % we only look for Y(t_{t+1})
            if problem == 0 || problem == 1
                Fk =  blkdiag(Rk_tr, Rk_tr)*fk;           
            elseif problem == 2
                Fk =  fk;
            elseif problem == 3
                Fk = zeros(6, 1);
            end
            Wm = (M + ht/2*K)*zm - (M - ht/2*K)*Yk +...
                    ht/4*(Q_yk*Yk + (Q_yk + Qdagger_yk)*zm + Q_zm*zm) +...
                    ht*Z*Fk;
            JWm = M + ht/2*K + ht/4*(Q_yk + Qdagger_yk) +...
                ht/4*(Q_zm + Qdagger_zm);
        end
        
        zxm1 = zxm - JWm\Wm; 

        %%% M.A. code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rel_err = (zxm1 - zxm)./zxm;                                  %
        nan_or_inf = find( isnan(rel_err) + isinf(rel_err) ...        %
            + (abs(zxm)<=tolzero) );                                  %
        rel_err(nan_or_inf) = 0;                                      %
        if (norm(rel_err,inf) <= reltolX) && (norm(Wm,inf) <= tolF)   %
            break                                                     %
        end                                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

        % update the variable zxm
        zxm = zxm1;
        zm = zxm1(1:Nf, 1);
        if approx_rot == 0 || approx_rot == 1
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
    
    % we just computed R(0, t_{k+1})
    if approx_rot == 0 || approx_rot == 1
        H(:, kk+1, 1) = xm;
        R_ini(:, :, kk+1) = func_quat2rotm(xm);

        % compute p(0, t_{k+1}):
        idx = 1:kk+1;
        v1_k1 = Y_phys(NNB(1:3, 1), kk+1);
        Bxt_ini(:, kk+1) = R_ini(:, :, kk+1)*v1_k1;
        p_ini(1, kk+1) = p_ini(1, 1) + trapz(t(idx), Bxt_ini(1, idx), 2);
        p_ini(2, kk+1) = p_ini(2, 1) + trapz(t(idx), Bxt_ini(2, idx), 2);
        p_ini(3, kk+1) = p_ini(3, 1) + trapz(t(idx), Bxt_ini(3, idx), 2);
    end

    %%% now we recover the position of the beam at this time k+1 %%%
    % Careful!! this scheme for computing p and R does not seem to converge
    % for problem == 0, 1, 2. Hence, latter on we will use another method
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
if plot_boundary == true
    if problem == 0 || problem == 1
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
end

%% Another way recover centerline
% Since solving for p and R using the TSolve (i.e. using velocities) does
% not seem to give a convergent scheme, here we use another method.
% We use XSolve (i.e. the strains) with p(0, t_k), R(0, t_k) as 'initial
% values'
centerline_mode = "XSolve";
[p2, R2] = recover_position(p(:, :, 1), R(:, :, :, 1), Y_phys, NNB, centerline_x_scheme, centerline_mode, x, t, flexMat, kap);
% [p2, R2] = recover_position(p_ini(:, :), R_ini(:, :, :), Y_phys, NNB, centerline_x_scheme, centerline_mode, x, t, flexMat, kap);
% [p2, R2] = recover_position(pD, RD, Y_phys, NNB, centerline_x_scheme, centerline_mode, x, t, flexMat, kap);

%% permute space and time for the variable p
p = permute(p, [1, 3, 2]); % for plotting we need to change the order of the space and time indexes

%% another way of computing p and R
% centerline_mode = "TSolve";      % or "XSolve"
% [p, R] = recover_position(p0, R0, Y_phys, NNB, centerline_t_scheme, centerline_mode, x, t, flexMat, kap);
%%% [p, R] = recover_position(p0_true, R0_true, y_true, NNB_true, centerline_t_scheme, centerline_mode, xNew, t, flexMat, kap);

%%% ----------- plot arclength ----------- %%%
centerline_mode = "TSolve";
file_name_arclength = ['fig/ARCLEN_p_', linNonlin, '_', centerline_mode, '.pdf'];
plot_arclen(p, x, t, centerline_mode, file_name_arclength);
centerline_mode = "XSolve";
file_name_arclength = ['fig/ARCLEN_p2_', linNonlin, '_', centerline_mode, '.pdf'];
plot_arclen(p2, x, t, centerline_mode, file_name_arclength);
%%% --------------------------------------- %%%

%%% ----------- plot centerline ----------- %%%
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
        f_c1 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        centerline_t = 1:(5*fact):Nt;  % time at which we plot the centerline
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        titleCenterline = 'Flying spaghetti problem 2D p2';
        f_c2 = plotCenterline(p2, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        
    elseif problem == 1
        titleCenterline = 'Flying spaghetti problem 3D';
        
%         centerline_t = fix([0, 2, 3, 3.8, 4.4, 5, 5.5, 5.8, 6.1, 6.5]*10*fact+1);
%         viewCent = [0, 90];
%         %viewCent = [90, 0];
%         locLegend = 'northeastoutside';
%         f_c1 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx);
% 
%         centerline_t = fix([0, 2.5, 3.5, 3.8, 4.5]*100+1);
%         viewCent = [90, 0];
%         locLegend = 'southwest';
%         f_c2 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx);
%         camroll(90);
        
        centerline_t = 1:(4*fact):Nt;
        viewCent = [0, 90]; % [0, 0] for x z plan
        locLegend = 'northeastoutside';
        f_c3 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
    
        titleCenterline = 'Flying spaghetti problem 3D p2';
        centerline_t = 1:(4*fact):Nt;
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c4 = plotCenterline(p2, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
    
    elseif problem == 2
        titleCenterline = 'Rotating arm p fin';
        centerline_t = (5*10*fact+1):(5*fact):(9*10*fact+1);
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c1 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        titleCenterline = 'Rotating arm p debut';
        centerline_t = 1:(5*fact):(5*10*fact+1);
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c2 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        
        titleCenterline = 'Rotating arm p2 fin';
        centerline_t = (5*10*fact+1):(5*fact):(9*10*fact+1);
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c3 = plotCenterline(p2, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        titleCenterline = 'Rotating arm p2 debut';
        centerline_t = 1:(5*fact):(5*10*fact+1);
        viewCent = [0, 90];
        locLegend = 'northeastoutside';
        f_c4 = plotCenterline(p2, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
    elseif problem == 3
        titleCenterline = 'Feedback control p - view 1';
        centerline_t = 1:(15*fact):Nt;
        locLegend = "southeast";
        viewCent = [96, 10];
        f_c1 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        titleCenterline = 'Feedback control p - view 2';
        locLegend = "northwest";
        viewCent = [-37.5, 30];
        f_c2 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        
        titleCenterline = 'Feedback control p2 - view 1';
        centerline_t = 1:(15*fact):Nt;
        locLegend = "southeast";
        viewCent = [96, 10];
        f_c3 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
        
        titleCenterline = 'Feedback control p2 - view 2';
        locLegend = "northwest";
        viewCent = [-37.5, 30];
        f_c4 = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t);
    end
end
%%% -------------------------------------- %%%
disp('End.')



%% definition of Q
function [res, resDagger] = Q(Y, P1, P2, P3, Pdagger1, Pdagger2, Pdagger3, Ni, Ne, Nf, NNB, diagonal, problem)
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
        if problem == 0 || problem == 1 || problem == 2
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
        elseif problem == 3
            for pp=1:6
                for ee = 2:Ne
                    res = res + Y(NNB(pp, 2*ee-1))*cell2mat(P1(pp, ee)) ...
                              + Y(NNB(pp, 2*ee))*cell2mat(P2(pp, ee)) ...
                              + Y(NNB(pp, 2*ee+1))*cell2mat(P3(pp, ee));
                    resDagger = resDagger + Y(NNB(pp, 2*ee-1))*cell2mat(Pdagger1(pp, ee)) ...
                                          + Y(NNB(pp, 2*ee))*cell2mat(Pdagger2(pp, ee)) ...
                                          + Y(NNB(pp, 2*ee+1))*cell2mat(Pdagger3(pp, ee));
                end
                ee=1;
                res = res + Y(NNB(pp, 2*ee))*cell2mat(P1(pp, ee)) ...
                          + Y(NNB(pp, 2*ee+1))*cell2mat(P2(pp, ee));
                resDagger = resDagger + Y(NNB(pp, 2*ee))*cell2mat(Pdagger1(pp, ee)) ...
                                      + Y(NNB(pp, 2*ee+1))*cell2mat(Pdagger2(pp, ee));
            end
            for pp=7:Ni
                for ee=1:Ne
                    res = res + Y(NNB(pp, 2*ee-1))*cell2mat(P1(pp, ee)) ...
                              + Y(NNB(pp, 2*ee))*cell2mat(P2(pp, ee)) ...
                              + Y(NNB(pp, 2*ee+1))*cell2mat(P3(pp, ee));
                    resDagger = resDagger + Y(NNB(pp, 2*ee-1))*cell2mat(Pdagger1(pp, ee)) ...
                                          + Y(NNB(pp, 2*ee))*cell2mat(Pdagger2(pp, ee)) ...
                                          + Y(NNB(pp, 2*ee+1))*cell2mat(Pdagger3(pp, ee));
                end
            end
        end
    end
end