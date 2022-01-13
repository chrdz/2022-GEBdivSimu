function plot1position(p, R, Nx)
% plot position of the beam at a specific time

%%% compute the cross sections:
r = 0.1; % radius
% border of the cross section is s \mapsto (0, r*cos(s), r*sin(s))
% for t in [0, 2*pi)
% position in space is thus: p0(x) + R0(x)*(0, r*cos(s), r*sin(s))
% for x in [0, \ell], s in [0, 2*pi)
Ns = 25;
s_list = linspace(0,2*pi,Ns);
circ = zeros(3, Ns, Nx);
for kk = 1:Nx
    for pp = 1:Ns
        s = s_list(pp);
        cross_pos = [0; r*cos(s); r*sin(s)];
        circ(:, pp, kk) = p(:, kk) + R(:, :, kk)*cross_pos;
    end
end

%%% plot of p0:
f = figure();
plot3(p(1, :), p(2, :), p(3, :), 'lineWidth', 2, 'Color','r');
hold on;
%%% plot cross sections:
for kk = 1:Nx
    axe_circ = plot3(circ(1, :, kk), circ(2, :, kk), circ(3, :, kk), 'linewidth', 2, 'Color','b');
    axe_circ.Color(4) = 0.5;
end
view(gca, 96,10);
axis equal;
grid on;
xlabel('X', 'Interpreter', 'latex', 'fontsize', 20);
ylabel('Y', 'Interpreter', 'latex', 'fontsize', 20);
zlabel('Z', 'Interpreter', 'latex', 'fontsize', 20);
title('Initial position of the beam', 'Interpreter', 'latex', 'fontsize', 20);
%exportgraphics(f,'fig/initial_centerline.pdf','ContentType','vector')


end

