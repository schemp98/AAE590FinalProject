function dxdt = linearized_rel_orbital_dynamics(t, x, mu, a, e, controller_type,B,R,Q)
% global u_history_global t_history_global
    % Evaluate orbital parameters at time t
    [rT, omegaT, omegaT_dot] = kepler_orbital_elements_eval(t, mu, a, e);
    [rTi, omegaTi, omegaT_doti] = kepler_orbital_elements_eval(0, mu, a, e);

    %Construct Dynamics (linearized)
    % A = zeros(6,6);
    % A(1:3, 4:6) = eye(3);
    % A(4,:) = [3*omegaT^2, omegaT_dot, 0, 0, 2*omegaT, 0];
    % A(5,:) = [-omegaT_dot, 0, 0, -2*omegaT, 0, 0];
    % A(6,:) = [0, 0, -omegaT^2, 0, 0, 0];

    rC = [rT + x(1); x(2); x(3)];
    normrC= norm(rC);
    mu_over_rC3 = -mu/(normrC^3);

    %Construct Dynamics (nonlinear)
    A = zeros(6,6);
    A(1:3, 4:6) = eye(3);
    A(4,:) = [ mu_over_rC3 + omegaT^2, omegaT_dot, 0, 0, 2*omegaT, 0];
    A(5,:) = [-omegaT_dot, mu_over_rC3 + omegaT^2, 0, -2*omegaT, 0, 0];
    A(6,:) = [0, 0, mu_over_rC3, 0, 0, 0];

    K = calculateControllerGainfunction(t, mu, a, e, controller_type,B,R,Q);
    
    u = K*x;

    dxdt = A*x + B*u;

    % %collect control history for plotting later
    % u_history_global(:,end+1) = u;
    % t_history_global = [t_history_global; t];
end


