function test_transfoIni(p0, R0, Y0, NNB, x, kap)
% test if we can recover p0, R0 from z0 via transfo 

Nx = length(x);
hx = x(2)-x(1);

q0_b = zeros(4, Nx);                    % quaternion
q0_b(:, 1) = rotm2quat(R0(:, :, 1));    % 'ini' data at x=0
R0_b = zeros(3, 3, Nx);                 % rotation
R0_b(:, :, 1) = quat2rotm(q0_b(:, 1)'); % 'ini' data at x=0

for kk=1:Nx-1     % mid point: find the quaternions
    z02_k = Y0(NNB(10:12, kk), 1);
    z02_k1 = Y0(NNB(10:12, kk+1), 1);
    Axt = func_U( (z02_k + z02_k1)/2 + kap ); 
    q0_b(:, kk+1) = ( eye(4) - hx/2 * Axt )\(( eye(4) + hx/2*Axt)*q0_b(:, kk));  
    R0_b(:, :, kk+1) = quat2rotm(q0_b(:, kk+1)'); % corresp. rotation
end

p0_b = zeros(3, Nx);    % position of centerline
p0_b(:, 1) = p0(:, 1);  % 'ini' data at x=0
MM = zeros(3, Nx); 
z01_1 = Y0(NNB(7:9, 1), 1);
MM(:, 1) = R0_b(:, :, 1)*(z01_1+[1; 0; 0]); 
for kk = 2:Nx      % integrate to find position of centerline   
    z01_k = Y0(NNB(7:9, kk), 1);
    MM(:, kk) = R0_b(:, :, kk)*(z01_k + [1; 0; 0]);
    idx = 1:kk;
    p0_b(1, kk) = p0_b(1, 1) + trapz(x(idx), MM(1, idx), 2);
    p0_b(2, kk) = p0_b(2, 1) + trapz(x(idx), MM(2, idx), 2);
    p0_b(3, kk) = p0_b(3, 1) + trapz(x(idx), MM(3, idx), 2);
end

%%% plot of p0 and p0_b:
%     figure();
%     plot3(p0(1, :), p0(2, :), p0(3, :), 'lineWidth', 2); hold on;
%     grid on;
%     plot3(p0_b(1, :), p0_b(2, :), p0_b(3, :), ':r', 'lineWidth', 2);

plot1position(p0_b, R0_b, Nx)

end

