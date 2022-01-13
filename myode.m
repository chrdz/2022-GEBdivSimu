function dydt = myode(tau,q, angular_velocity, t_list)
% function for ode 45 to in the scheme to recover the quaternions and then
% the position of the centerline
    f = interp1(t_list,transpose(angular_velocity),tau); % Interpolate the data set (ft,f) at time t
    dydt = func_U(transpose(f))*q; % Evaluate ODE at time t
end