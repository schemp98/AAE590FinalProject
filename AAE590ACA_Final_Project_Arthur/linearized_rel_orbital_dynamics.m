function dxdt = linearized_rel_orbital_dynamics(t, x, mu, a, e)
    % Evaluate orbital parameters at time t
    [omegaT, omegaT_dot] = kepler_orbital_elements_eval(t, mu, a, e);

    %Construct Dynamics
    A = zeros(6,6);
    A(1:3, 4:6) = eye(3);
    A(4,:) = [3*omegaT^2, omegaT_dot, 0, 0, 2*omegaT, 0];
    A(5,:) = [-omegaT_dot, 0, 0, -2*omegaT, 0, 0];
    A(6,:) = [0, 0, -omegaT^2, 0, 0, 0];

    %implement dynamics (debugging with CWH)
    % n = sqrt(mu/a^3);
    % A = [zeros([3 3]) eye(3); 3*n^2 0 0 0 2*n 0;...
    %     zeros([1 3]) -2*n 0 0; 0 0 -n^2 0 0 0];
    % B = [zeros([3 3]); eye(3)];

    B = [zeros(3,3); eye(3)];

    
    u = zeros(3,1); %Eventually replace with SDRE input directly here? 

    dxdt = A*x + B*u;
end
