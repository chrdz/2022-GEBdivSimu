function [p, R] = recover_position(p0, R0, Y, NNB, centerline_scheme, type_centerline, x, t, flexMat, kap)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    % arguments: 
    % p0, R0 initial position of the centerline
    % Y solution 
    % NNB node numbers 
    % ht time step
    % centerline_t_scheme which scheme we will use
    % returns:
    % p, R position of the beam through time
    
    Nx = length(x);
    Nt = length(t);
    hx = x(2) - x(1);
    ht = t(2) - t(1);

    %% recover position of the centerline using velocities
    if type_centerline == 'TSolve'
        %%% WITH TIME : using the first 6 eq. of the transformation %%%
        disp('Recovering the position of the centerline (using the velocities)..')
        tic

        H = zeros(4, Nt, Nx);     % quaternions
        R = zeros(3, 3, Nt, Nx);  % rotation matrices
        p = zeros(3, Nt, Nx);     % positions
        for kk = 1:Nx                                 % initial data at t=0
            H(:, 1, kk) = rotm2quat(R0(:, :, kk))';   % quat angle at t=0
            R(:, :, 1, kk) = R0(:, :, kk);            % angle at t=0
            p(:, 1, kk) = p0(:, kk);                  % position at t=0
        end
        % we have y(x,t) containing the velocities and strains
        % here we will just use the velocities
        % the velocoties are the first six components of y
        velocities = zeros(6, Nt, Nx); 
        for kk = 1:Nx               % for all x
            for nn = 1:Nt%-1         % for all t
                for ii = 1:6
                    velocities(ii, nn, kk) = Y(NNB(ii, kk), nn); % extract velocities
                end
            end
        end

        if centerline_scheme == 0       %%% mid point rule
            for kk = 1:Nx                 % for all x
                for nn = 1:Nt-1           % for all t
                    Axt = func_U( (velocities(4:6, nn, kk)+velocities(4:6, nn+1, kk) )/2);  % uses angular velocities  
                    H(:, nn+1, kk) = (eye(4)-ht/2*Axt)\((eye(4) + ht/2*Axt)*H(:,nn, kk)); % GOOD
                    %disp(norm(H(:, nn+1, kk)))  % test
                    %R(:, :, nn+1, kk) = quat2rotm(H(:, nn+1, kk)');
                    R(:, :, nn+1, kk) = func_quat2rotm(H(:, nn+1, kk));
                end
            end
        elseif centerline_scheme == 1   %%% explicit Euler (does not work properly)
            for kk = 1:Nx                 % for all x
                for nn = 1:Nt-1           % for all t
                    Axt = func_U( velocities(4:6, nn, kk));  % uses angular velocities  
                    H(:, nn+1, kk) = (eye(4) + ht*Axt)*H(:, nn, kk);
                    H(:, nn+1, kk) = H(:, nn+1, kk)/norm(H(:, nn+1, kk));
                    %R(:, :, nn+1, kk) = quat2rotm(H(:, nn+1, kk)');
                    R(:, :, nn+1, kk) = func_quat2rotm(H(:, nn+1, kk));
                end
            end
        elseif centerline_scheme == 2   %%% ode 45
            tspan = t;
            opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
        for kk = 2:Nx
            ic = H(:, 1, kk); 
            [t_ode,y_ode] = ode45(@(tau,q) myode(tau,q,velocities(4:6, :, kk), tspan), tspan, ic, opts);
            H(:, :, kk) = transpose(y_ode);
            for nn = 1:Nt
                %R(:, :, nn, kk) = quat2rotm(H(:, nn, kk)');
                R(:, :, nn, kk) = func_quat2rotm(H(:, nn, kk));
            end
        end
        elseif centerline_scheme == 3   %%% implicit euler
            for kk = 1:Nx                 % for all x
                for nn = 1:Nt-1           % for all t
                    Axt = func_U( velocities(4:6, nn+1, kk));  % uses angular velocities  
                    %abs(eig((eye(4) - ht*Axt)\eye(4)))   % test
                    H(:, nn+1, kk) = (eye(4) - ht*Axt)\H(:, nn, kk);
                    %R(:, :, nn+1, kk) = quat2rotm(H(:, nn+1, kk)');
                    R(:, :, nn+1, kk) = func_quat2rotm(H(:, nn+1, kk));
                end
            end
        else                              %%% Z. paper (does not work properly)
            for kk = 1:Nx                 % for all x
                for nn = 1:Nt-1           % for all t
                    Wm = (velocities(4:6, nn, kk)+velocities(4:6, nn+1, kk) )/2;
                    AAA = (func_U_right(Wm) + eye(4)*4/ht)*4*ht^2/(16+ht^2*norm(Wm)^2);            
                    qqq0 = H(1, nn, kk);
                    qqq = H(2:4, nn, kk);
                    bbb = [1/ht*qqq0 - 1/4*transpose(qqq)*Wm;...
                        1/ht*qqq + 1/4*(qqq0*Wm + hat(Wm)*qqq)];
                    H(:, nn+1, kk) = AAA*bbb/norm(AAA*bbb);
                    %disp(norm(H(:, nn+1, kk)))   % test
                    R(:, :, nn+1, kk) = func_quat2rotm(H(:, nn+1, kk));
                end
            end    
        end

        Bxt = zeros(3, Nt, Nx); 
        for kk = 1:Nx
            Bxt(:, 1, kk) = R(:, :, 1, kk)*(velocities(1:3, 1, kk));
            for nn = 2:Nt
                idx = 1:nn;
                Bxt(:, nn, kk) = R(:, :, nn, kk)*(velocities(1:3, nn, kk));
                p(1, nn, kk) = p(1, 1, kk) + trapz(t(idx), Bxt(1, idx, kk), 2);
                p(2, nn, kk) = p(2, 1, kk) + trapz(t(idx), Bxt(2, idx, kk), 2);
                p(3, nn, kk) = p(3, 1, kk) + trapz(t(idx), Bxt(3, idx, kk), 2);
            end
        end
        p = permute(p, [1, 3, 2]); % for plotting we need to change the order of the space and time indexes
        R = permute(R, [1, 2, 4, 3]);
        
        toc
    else
        %% recover position of the centerline using internal forces and moments
        disp('Recovering the position of the centerline (using the internal forces and moments / strains)..')
        tic

        %%% WITH SPACE: using the last 6 eq. of the transformation  %%%
        H = zeros(4, Nx, Nt);     % quaternions
        R = zeros(3, 3, Nx, Nt);  % rotation matrices
        p = zeros(3, Nx, Nt);     % positions
        for nn = 1:Nt                                 % clamped at all times
            H(:, 1, nn) = rotm2quat(R0(:, :, nn))';    % quat angle at x=0
            R(:, :, 1, nn) = R0(:, :, nn);             % angle at x=0
            p(:, 1, nn) = p0(:, nn);                   % position at x=0
        end
        % we have y(x,t) containing the velocities and strains
        % here we will just use the strains
        % the strains are the last six components of y left-multiplied by flexMat
        strains = zeros(6, Nx, Nt); 
        for nn = 1:Nt            % for all t
            for kk = 1:Nx%-1      % for all x
                forces = zeros(6, 1);  
                for ii = 7:12
                    forces(ii-6, 1) = Y(NNB(ii, kk), nn);    % extract the forces
                end
                strains(:, kk, nn) = flexMat*forces;         % left-mult by flexMat
            end
        end

        if centerline_scheme == 0  %%% mid point rule
            for nn = 1:Nt            % for all t
                for kk = 1:Nx-1      % for all x
                    Axt = func_U(strains(4:6, kk, nn)/2 + strains(4:6, kk+1, nn)/2 + kap );
                    H(:, kk+1, nn) = (eye(4) - hx/2*Axt)\((eye(4) + hx/2*Axt)*H(:, kk, nn));
                    %R(:, :, kk+1, nn) = quat2rotm(H(:, kk+1, nn)');
                    R(:, :, kk+1, nn) = func_quat2rotm(H(:, kk+1, nn));
                end
            end
        else                         %%% euler explicit
            for nn = 1:Nt            % for all t
                for kk = 1:Nx-1      % for all x
                    Axt = func_U( strains(4:6, kk, nn) + kap );  
                    H(:, kk+1, nn) = (eye(4) + hx*Axt)*H(:, kk, nn);
                    %R(:, :, kk+1, nn) = quat2rotm(H(:, kk+1, nn)');
                    R(:, :, kk+1, nn) = func_quat2rotm(H(:, kk+1, nn));
                end
            end
        end

        Bxt = zeros(3, Nx, Nt); 
        for nn = 1:Nt
            Bxt(:, 1, nn) = R(:, :, 1, nn)*(strains(1:3, 1, nn) + [1; 0; 0]);
            for kk = 2:Nx
                idx = 1:kk;
                Bxt(:, kk, nn) = R(:, :, kk, nn)*(strains(1:3, kk, nn) + [1; 0; 0]);
                p(1, kk, nn) = p(1, 1, nn) + trapz(x(idx), Bxt(1, idx, nn), 2);
                p(2, kk, nn) = p(2, 1, nn) + trapz(x(idx), Bxt(2, idx, nn), 2);
                p(3, kk, nn) = p(3, 1, nn) + trapz(x(idx), Bxt(3, idx, nn), 2);
            end
        end

        toc
    end
        
end

