function rel_LVLH_hist = inertial_diff_to_LVLH(rT_hist, vT_hist, rC_hist)
% inertial_diff_to_LVLH
% Transforms the inertial difference vector between chaser and target
% into the target's LVLH frame.
%
% Inputs:
%   rT_hist - Nx3 matrix of target inertial positions [m]
%   vT_hist - Nx3 matrix of target inertial velocities [m/s]
%   rC_hist - Nx3 matrix of chaser inertial positions [m]
%
% Output:
%   rel_LVLH_hist - Nx3 matrix of relative position in LVLH frame [m]

    N = size(rT_hist,1);
    rel_LVLH_hist = zeros(N,3);

    for k = 1:N
        rT = rT_hist(k,:)';
        vT = vT_hist(k,:)';
        rC = rC_hist(k,:)';

        % Compute LVLH frame basis vectors
        rHat = rT / norm(rT);
        hVec = cross(rT, vT);
        hHat = hVec / norm(hVec);
        y_LVLH = -hHat;                % negative angular momentum direction
        z_LVLH = -rHat;                % negative radial direction
        x_LVLH = cross(y_LVLH, z_LVLH); % completes right-handed triad

        R_LVLH_to_I = [x_LVLH, y_LVLH, z_LVLH];

        % Transform inertial difference into LVLH frame
        rho_inertial = rC - rT;
        rho_LVLH = R_LVLH_to_I' * rho_inertial;  % transpose = inertial-to-LVLH rotation
        rel_LVLH_hist(k,:) = rho_LVLH';
    end
end
