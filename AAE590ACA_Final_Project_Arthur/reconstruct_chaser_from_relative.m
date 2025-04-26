function [rC_hist, vC_hist] = reconstruct_chaser_from_relative(t, x_rel_hist, rT_hist, vT_hist)

%candidate function for transforming relative state to inertial state

    N = length(t);
    rC_hist = zeros(N, 3);
    vC_hist = zeros(N, 3);

    for k = 1:N
        % Extract target state at time t(k)
        rT = rT_hist(k, :)';
        vT = vT_hist(k, :)';

        % Relative state in LVLH
        rho = x_rel_hist(k, 1:3)';
        rho_dot = x_rel_hist(k, 4:6)';

        % Construct LVLH basis
        rHat = rT / norm(rT);
        hVec = cross(rT, vT);
        hHat = hVec / norm(hVec);
        y_LVLH = -hHat;
        z_LVLH = -rHat;
       % x_LVLH = cross(y_LVLH, z_LVLH);
       x_LVLH = cross(z_LVLH, y_LVLH);

        % Rotation matrix: LVLH -> inertial
        R_LVLH_to_I = [x_LVLH, y_LVLH, z_LVLH];

        % Angular velocity of LVLH frame
        omega_LVLH = cross(rT, vT) / norm(rT)^2;

        % Transform position and velocity
        rho_I = R_LVLH_to_I * rho;
        rho_dot_I = R_LVLH_to_I * (rho_dot + cross(omega_LVLH, rho));

        % Reconstruct chaser state in inertial frame
        rC_hist(k, :) = (rT + rho_I)';
        vC_hist(k, :) = (vT + rho_dot_I)';
    end
end
