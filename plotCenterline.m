function f_c = plotCenterline(p, centerline_t, viewCent, locLegend, titleCenterline, Nx, t)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


    f_c = figure();
    first = centerline_t(1);
    h = plot3(p(1, :, first), p(2, :, first), p(3, :, first), 'g', 'lineWidth', 2); 
    hold on;
    %h.Color(4) = 0.5;
    grid on;
    sec = centerline_t(2);
    h = plot3(p(1, :, sec), p(2, :, sec), p(3, :, sec), 'k', 'lineWidth', 1);
    for nn = centerline_t(3:end-1)
        h = plot3(p(1, :, nn), p(2, :, nn), p(3, :, nn), 'k', 'lineWidth', 1);
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    last = centerline_t(end);
    plot3(p(1, :, last), p(2, :, last), p(3, :, last), 'r', 'lineWidth', 2);
    %plot3(pf(1, :), pf(2, :), pf(3, :), '--g', 'lineWidth', 3);

    pp = permute(p, [1, 3, 2]);
    plot3(pp(1, first:last, 1), pp(2, first:last, 1), pp(3, first:last, 1), '--b');
    plot3(pp(1, first:last, Nx), pp(2, first:last, Nx), pp(3, first:last, Nx), '--r');

    axis equal;
    view(gca, viewCent);
    legend(['$t = ', num2str(t(first)), '$'], '$t$', ['$t = ', num2str(t(last)), '$'],...
        '$x=0$', '$x=\ell$','Interpreter','latex',...
        'Location', locLegend, 'fontsize', 12);
    xlabel('X', 'Interpreter', 'latex');
    ylabel('Y', 'Interpreter', 'latex');
    zlabel('Z', 'Interpreter', 'latex');
    title(titleCenterline, 'Interpreter', 'latex', 'fontsize', 12);
    %exportgraphics(f_c,'flying_spaghetti.pdf','ContentType','vector')


end

