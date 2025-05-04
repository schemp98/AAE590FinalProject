function K = calculateControllerGainfunction(t, x, mu, a, e, controller_type,B,R,Q)

[rT, omegaT, omegaT_dot] = kepler_orbital_elements_eval(t, mu, a, e);


    %Construct Dynamics
    % A = zeros(6,6);
    % A(1:3, 4:6) = eye(3);
    % A(4,:) = [3*omegaT^2, omegaT_dot, 0, 0, 2*omegaT, 0];
    % A(5,:) = [-omegaT_dot, 0, 0, -2*omegaT, 0, 0];
    % A(6,:) = [0, 0, -omegaT^2, 0, 0, 0];

    rC = [rT + x(1); x(2); x(3)];
    normrC= norm(rC);
    mu_over_rC3 = mu/(normrC^3);

    % A = zeros(6,6);
    % A(1:3, 4:6) = eye(3);
    % A(4,:) = [ mu_over_rC3 + 2*omegaT^2, omegaT_dot, 0, 0, 2*omegaT, 0];
    % A(5,:) = [-omegaT_dot, mu_over_rC3 + omegaT^2, 0, -2*omegaT, 0, 0];
    % A(6,:) = [0, 0, mu_over_rC3, 0, 0, 0];

        %Construct Dynamics (in class linearized EOMs)
    A = zeros(6,6);
    A(1:3, 4:6) = eye(3);
    A(4,:) = [ 2*mu_over_rC3 + omegaT^2, omegaT_dot, 0, 0, 2*omegaT, 0];
    A(5,:) = [-omegaT_dot, omegaT^2 - mu_over_rC3, 0, -2*omegaT, 0, 0];
    A(6,:) = [0, 0, -mu_over_rC3, 0, 0, 0];

    switch controller_type
        case 'SDRE'
            %implement SDRE
            %define hamiltonian matrix
            HM = [A,    -B*inv(R)*B';
                  -Q,        -A'    ];
        
            % [V,D] = eig(HM);
            % eigvals = real(diag(D));
            % 
            % n = 0;
            % for i = 1:length(eigvals)
            %     if eigvals(i) > 0
            %         n = n+1;
            %     end
            % end
            % 
            % if n ~= 6
            %     error('wrong number of eigvals with negative real part')
            % end
        
            [Y,X] = ric_schr(HM);
    
            P = X*inv(Y);
        
            K = -inv(R)* B'* P;
        
            % riccati_val = P*A + A'*P - P*B*inv(R)*B'*P + Q;
        
            % u =  K * x;

        case 'LQR'
            %use the CWH definition of A to make it LTI? 
            n = sqrt(mu/a^3);
           A_lqr = [zeros([3 3]) eye(3); 3*n^2 0 0 0 2*n 0;...
         zeros([1 3]) -2*n 0 0; 0 0 -n^2 0 0 0];
         [K,~,~] = lqr(A_lqr,B,Q,R);

         K = -K;
            % u = -K_lqr * x;

        otherwise
            error('Unknown controller type: %s', controller_type);
    end