
if problem == 0
    movie_name = 'flying_spaghetti_20_2D';
    viewCent = [0, 90];
    locLegend = 'northeastoutside';
    title_pic = 'flying_spaghetti_20_2D.pdf';
elseif problem == 1
    movie_name = 'flying_spaghetti_20_3D';
    locLegend = 'northeastoutside';
    title_pic = 'flying_spaghetti_20_3D.pdf';
elseif problem == 2
    titleCenterline = 'Rotating arm p2';
    movie_name = 'rotating_arm_15_2D';
    viewCent = [0, 90];
    locLegend = 'northeastoutside';
    title_pic = 'rotating_arm_15_2D.pdf';
end

movie_name = 'test';

fig = movieCenterline(p2, R2, 1:numel(t), viewCent, locLegend, titleCenterline, t, movie_name);
title('Random title', 'Color','none')
%exportgraphics(fig,title_pic,'ContentType','vector');



%ffmpeg -i rotating_arm_15_2D.mp4 rotating_arm_15_2D.gif
%sips -s format png flying_spaghetti_20_2D.pdf --out flying_spaghetti_20_2D.png