function res = mom1_book(t)
% function M_1
% to define the external force

    if t < 5
        res = 1.5/5;
    else
        res = 0;
    end
    

end


% function res = mom1_book(t)
% % function M_1
% % to define the external force
% 
%     tmax = 5;
%     x0 = 0.5*tmax;   % center
%     a = 0.2*tmax;    % width
%     c0 = 1.5/5;       % magnitude
%     f = @(tt) c0 * exp( -1/a^2*((tt-x0).^2) ) ; % bump function
%     res = f(t);
%    
% end