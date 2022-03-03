% function res = f_ext(t)
% %UNTITLED5 Summary of this function goes here
% %   Detailed explanation goes here
% res = zeros(6, 1);
% res(1, 1) = for1(t);
% res(4, 1) = mom1(t);
% res(5, 1) = mom2(t);
% end


function res = f_ext(t, problem)
%planar case
    res = zeros(6, 1);

    if problem == 0
        res(1, 1) = -for1_planar(t);
        res(6, 1) = mom1_planar(t);
    elseif problem == 1
        res(1, 1) = -for1(t);
        %simo version:
        res(5, 1) = -mom2(t); % adding minus sign does note seem to change anything
        res(6, 1) = mom1(t);
        
%         %hesse verison TEST:
%         res(5, 1) = -mom1(t);
%         res(6, 1) = -mom2(t);
        
    elseif problem == 2
        res(6, 1) = mom1_book(t);
    end
end

% function res = f_ext(t)
% %UNTITLED5 Summary of this function goes here
% %   Detailed explanation goes here
% res = zeros(6, 1);
% res(1, 1) = -0.6*for1(t);
% res(2, 1) = 0.8*for1(t);
% % res(5, 1) = 0.8*mom1(t);
% % res(6, 1) = -0.6*mom1(t)+ mom2(t);
% 
% res(5, 1) = 0.8*mom1(t);
% res(6, 1) = -0.6*mom1(t)+ mom2(t);
% end

% function res = f_ext(t)
% %planar case
% res = zeros(6, 1);
% res(1, 1) = -0.6*for1_planar(t);
% res(2, 1) = 0.8*for1_planar(t);
% %res(1, 1) = for1_planar(t);
% 
% res(4, 1) = mom1_planar(t);
% 
% % res(6, 1) = mom1_planar(t);
% % res(5, 1) = mom1_planar(t);
% end


%% DNA
% % % % % function res = f_ext(t)
% % % % % %UNTITLED5 Summary of this function goes here
% % % % % %   Detailed explanation goes here
% % % % % res = zeros(6, 1);
% % % % % res(1, 1) = -0.6*for1(t);
% % % % % res(2, 1) = 0.8*for1(t);
% % % % % res(5, 1) = 0.8*mom1(t);
% % % % % res(6, 1) = -0.6*mom1(t) + mom2(t);
% % % % % end


