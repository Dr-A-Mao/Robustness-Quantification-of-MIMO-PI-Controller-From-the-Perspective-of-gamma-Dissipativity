function [J_F0,J_G0] = get_Jacobians(phi_c,theta_c,J_x,J_y,J_z,type)
    if nargin < 1
        phi_c =  pi/3;
        theta_c = pi/6;

        J_x = 9496  * 1.3558;
        J_y = 55814 * 1.3558;
        J_z = 63100 * 1.3558;

        type = 2;
    end
    % UAV attitude control
    if type == 1
        J_F0 = zeros(6);
        J_F0(1:3,4:6) = [1,tan(theta_c) * sin(phi_c),tan(theta_c) * cos(phi_c);...
                         0,cos(phi_c),-sin(phi_c);...
                         0,sin(phi_c)/cos(theta_c),cos(phi_c)/cos(theta_c)];
        J_G0 = [zeros(3);diag([-1/J_x,-1/J_y,-1/J_z])];
    end
    % Fixed-wing UAV guidance law
    if type == 2
        g = 9.81; V = 25; gamma_c = pi/12;
        J_F0 = [0,0;0,g/V * sin(gamma_c)];
        J_G0 = diag([-g/V,-g/V]);
    end
end