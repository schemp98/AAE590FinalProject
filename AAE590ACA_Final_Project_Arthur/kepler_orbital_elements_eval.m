function [rT, omegaT, omegaT_dot] = kepler_orbital_elements_eval(t, mu, a, e)
% global drdt_log t_log omegaT_log omegaT_dot_log
    % Constants
    n = sqrt(mu / a^3); %mean motion
    M = n * t; %mean anomaly

    %Get Eccentric anomaly
    E = mean_to_ecc_anomaly_2(M,e);

    %true anom
    theta = 2 * atan2(sqrt(1 + e)*sin(E/2), sqrt(1 - e)*cos(E/2));

    p = a * (1 - e^2);

    rT = p / (1 + e * cos(theta));

    h = sqrt(mu * a * (1 - e^2));
    % p_test = h^2/mu;

    dr_dtheta = (p * e * sin(theta)) / (1 + e * cos(theta))^2;

    h = sqrt(mu * p);
    dtheta_dt = h / rT^2;

    % z = sqrt(1-e^2);
    % 
    % dtheta_dt_test = z*(1 + e*cos(theta))^2/(z^3); %----debugging

    %dtheta_dt_test - dtheta_dt
    
    omegaT = sqrt(mu / rT^3);
    
    %from vallado for noncircular orbits
    %omegaT = sqrt(mu*p / (rT^4));

    omegaT_dot = -(3/2) * rT^(-5/2) * sqrt(mu) * dr_dtheta * dtheta_dt;

    %collecting variables for noncircular orbit diagnostics---------------
    %commenting out as its no longer needed

   % drdt_log(end+1) = dr_dtheta * dtheta_dt;
   % t_log(end+1) = t;
   % omegaT_log(end +1) = omegaT;
   % omegaT_dot_log(end +1) = omegaT_dot;

end
