function f_c = movieCenterline(p, R, t_idx, viewCent, locLegend, titleCenterline, t)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    Nt = numel(t);
    Nx = size(p, 2);
    ht = t(2) - t(1);
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
    Ns = 100;
    s_list = linspace(0,2*pi,Ns); % the angles 
    circ = zeros(3, Nx, Ns, Nt);
    for pp = 1:Ns
        s = s_list(pp); % one angle
        cross_pos = [0; r*cos(s); r*sin(s)];
        for nn = 1:Nt
            for kk = 1:Nx
                circ(:, kk, pp, nn) = p(:, kk, nn) + R(:, :, kk, nn)*cross_pos;
            end
        end
    end
    
    for nn = t_idx(1):(1/ht/10):t_idx(end)
        
        figure(2) %%% ADD THIS HERE otherwise it will plot in fig 1 for kk>1 (because of the hold off) 
        
        %%% plot cross sections:
        for pp = 1:Ns
            axe_circ = plot3(circ(1, :, pp, nn), circ(2, :, pp, nn),...
                circ(3, :, pp, nn), 'linewidth', 2, 'Color','b');
            axe_circ.Color(4) = 0.01;
            hold on;
        end
        plot3(p(1, :, nn), p(2, :, nn), p(3, :, nn), 'r', 'lineWidth', 1);
        for ss = 1:(1/ht/10):nn-1
            hist = plot3(p(1, :, ss), p(2, :, ss), p(3, :, ss), 'k:', 'lineWidth', 1);
            hist.Color(4) = 0.25;
        end
        
        pp = permute(p, [1, 3, 2]);
        plot3(pp(1, 1:nn, 1), pp(2, 1:nn, 1), pp(3, 1:nn, 1), ':b', 'lineWidth', 1.5);
        plot3(pp(1, 1:nn, Nx), pp(2, 1:nn, Nx), pp(3, 1:nn, Nx), ':r', 'lineWidth', 1.5);
        
        grid on;
        axis equal;
        %view(gca, viewCent);
        xlim([x_min, x_max]);
        ylim([y_min, y_max]);
        zlim([z_min, z_max]);
        xlabel('X', 'Interpreter', 'latex');
        ylabel('Y', 'Interpreter', 'latex');
        zlabel('Z', 'Interpreter', 'latex');
        title(titleCenterline, 'Interpreter', 'latex', 'fontsize', 12);
        hold off;

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

