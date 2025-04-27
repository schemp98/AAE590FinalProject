function dxdt = linearized_rel_orbital_dynamics(t, x, mu, a, e, controller_type)
global u_history_global t_history_global
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

    B = [zeros(3,3); eye(3)];

    %define Q and R
    Q = zeros(6,6);

    Q(1:3,1:3) = eye(3).*0.01;
    Q(4:6,4:6) = eye(3).*0.001;

    R = eye(3).*10^7;

    switch controller_type
        case 'SDRE'
            %implement SDRE
            %define hamiltonian matrix
            HM = [A,    -B*inv(R)*B';
                  -Q,        -A'    ];
        
            [V,D] = eig(HM);
            eigvals = real(diag(D));

            n = 0;
            for i = 1:length(eigvals)
                if eigvals(i) > 0
                    n = n+1;
                end
            end
        
            if n ~= 6
                error('wrong number of eigvals with negative real part')
            end
        
            [Y,X] = ric_schr(HM);
    
            P = X*inv(Y);
        
            k = -inv(R)* B'* P;
        
            riccati_val = P*A + A'*P - P*B*inv(R)*B'*P + Q;
        
            u =  k * x;

        case 'LQR'
            %use the CWH definition of A to make it LTI? 
            n = sqrt(mu/a^3);
           A_lqr = [zeros([3 3]) eye(3); 3*n^2 0 0 0 2*n 0;...
         zeros([1 3]) -2*n 0 0; 0 0 -n^2 0 0 0];
         [K_lqr,~,~] = lqr(A_lqr,B,Q,R);
        % 
        % A_lqr = zeros(6,6);
        % A_lqr(1:3, 4:6) = eye(3);
        % A_lqr(4,:) = [3*omegaTi^2, omegaT_doti, 0, 0, 2*omegaTi, 0];
        % A_lqr(5,:) = [-omegaT_doti, 0, 0, -2*omegaTi, 0, 0];
        % A_lqr(6,:) = [0, 0, -omegaTi^2, 0, 0, 0];
        % [K_lqr,~,~] = lqr(A_lqr,B,Q,R);

            u = -K_lqr * x;

        otherwise
            error('Unknown controller type: %s', controller_type);
    end

    dxdt = A*x + B*u;

    %collect control history for plotting later
    u_history_global(:,end+1) = u;
    t_history_global = [t_history_global; t];
end


