function f_c = movieCenterline(p, R, centerline_t, viewCent, locLegend, titleCenterline, Nx, t)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    Nt = numel(t);
    create_movie = true;
    
    if create_movie
        movie_name = 'test_movie';  % file name
%         vidfile = VideoWriter(movie_name,'Uncompressed AVI');
        vidfile = VideoWriter(movie_name,'MPEG-4');
        vidfile.FrameRate = 10;      % change this number to slow down or speed up the movie
        open(vidfile);
        fig=figure(2);
        set(fig,'color','w');
    end
    
    x_min = min(min(p(1, :, :)));
    x_max = max(max(p(1, :, :)));
    y_min = min(min(p(2, :, :)));
    y_max = max(max(p(2, :, :)));
    z_min = min(min(p(3, :, :)));
    z_max = max(max(p(3, :, :)));
    
    
%     r = 0.025; % radius
    r = 0.1; % radius
    % border of the cross section is s \mapsto (0, r*cos(s), r*sin(s))
    % for t in [0, 2*pi)
    % position in space is thus: p0(x) + R0(x)*(0, r*cos(s), r*sin(s))
    % for x in [0, \ell], s in [0, 2*pi)
    Ns = 25;
    s_list = linspace(0,2*pi,Ns); % the angles 
    circ = zeros(3, Ns, Nx, Nt);
    for pp = 1:Ns
        s = s_list(pp); % one angle
        cross_pos = [0; r*cos(s); r*sin(s)];
        for nn = 1:Nt
            for kk = 1:Nx
                circ(:, pp, kk, nn) = p(:, kk, nn) + R(:, :, kk, nn)*cross_pos;
            end
        end
    end
    
    for nn = centerline_t(1):10:centerline_t(end)
        
        figure(2) %%% ADD THIS HERE otherwise it will plot in fig 1 for kk>1 (because of the hold off) 
        
        plot3(p(1, :, nn), p(2, :, nn), p(3, :, nn), 'k', 'lineWidth', 1);
        hold on;
%         %%% plot cross sections:
%         for kk = 1:Nx
%             axe_circ = plot3(circ(1, :, kk, nn), circ(2, :, kk, nn),...
%                 circ(3, :, kk, nn), 'linewidth', 2, 'Color','b');
%             axe_circ.Color(4) = 0.1;
%         end
        
        
        grid on;
        axis equal;
        %view(gca, viewCent);
%         xlim([0, 50]);
%         ylim([0, 5]);
%         zlim([-5, 5]);
        xlim([x_min, x_max]);
        ylim([y_min, y_max]);
        zlim([z_min, z_max]);
        xlabel('X', 'Interpreter', 'latex');
        ylabel('Y', 'Interpreter', 'latex');
        zlabel('Z', 'Interpreter', 'latex');
        title(titleCenterline, 'Interpreter', 'latex', 'fontsize', 12);
        hold off;
%         view(2)
%         title(['t = ', num2str(time(kk)), ' [s]'])
%         xlabel 'x [m]'
%         ylabel 'y [m]'
%         axis tight
%         hcb = colorbar();
%         ylabel(hcb,'concentration of Cl [mg/l]')
%         colormap default
%         shading interp
%         axis equal
%         hold off

        if create_movie
            frame = getframe(fig);
            writeVideo(vidfile, frame);
        end

        pause(0.1);
    end 
    
    if create_movie
        close(vidfile);
    end

end

