function [r_eci, v_eci] = keplerian2cartesian(a, e, i, RAAN, argp, theta, mu)

%test function

% keplerian2cartesian - Convert Keplerian elements to ECI Cartesian state
%
% Inputs:
%   a     - semi-major axis [m]
%   e     - eccentricity
%   i     - inclination [deg]
%   RAAN  - right ascension of ascending node [deg]
%   argp  - argument of perigee [deg]
%   theta - true anomaly [deg]
%   mu    - gravitational parameter [m^3/s^2]
%
% Outputs:
%   r_eci - position vector in ECI [m]
%   v_eci - velocity vector in ECI [m/s]

    % Convert angles to radians
    i = deg2rad(i);
    RAAN = deg2rad(RAAN);
    argp = deg2rad(argp);
    theta = deg2rad(theta);

    % Compute orbital radius
    r = a * (1 - e^2) / (1 + e * cos(theta));

    % Position in perifocal frame
    r_pf = r * [cos(theta); sin(theta); 0];

    % Velocity in perifocal frame
    h = sqrt(mu * a * (1 - e^2));
    v_pf = (mu / h) * [-sin(theta); e + cos(theta); 0];

    % Rotation matrix from perifocal to ECI
    Rz_RAAN = [cos(RAAN), -sin(RAAN), 0;
               sin(RAAN),  cos(RAAN), 0;
               0,          0,         1];

    Rx_i = [1, 0, 0;
            0, cos(i), -sin(i);
            0, sin(i),  cos(i)];

    Rz_argp = [cos(argp), -sin(argp), 0;
               sin(argp),  cos(argp), 0;
               0,          0,         1];

    Q_pX = Rz_RAAN * Rx_i * Rz_argp;

    % Transform to ECI
    r_eci = Q_pX * r_pf;
    v_eci = Q_pX * v_pf;
end
