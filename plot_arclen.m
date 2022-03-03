function plot_arclen(p, x, t, centerline_mode, file_name_arclength)
% plot arclength

if centerline_mode == "XSolve"
    title_arclength = 'Arclength: by space integration using $z$';
elseif centerline_mode == "TSolve"
    title_arclength = 'ArcLength: by time integration using $v$';
end

Nt = length(t);
hx = x(2)-x(1);

arcLength_time = zeros(Nt, 1);
for nn = 1:Nt
    arcLength_time(nn, 1) = trapz(x, sqrt(sum((gradient(p(:, :, nn), hx)).^2))); % arclength
end
f_arc = figure();
plot(t, arcLength_time, 'lineWidth', 2);
xlabel('$t$','Interpreter','latex');
grid on
title(title_arclength, 'Interpreter', 'latex');
%print(f_arc, file_name_arclength,'-dpdf')
%exportgraphics(f_arc,file_name_arclength,'ContentType','vector');

end

