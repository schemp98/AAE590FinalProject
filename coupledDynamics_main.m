

addpath("AAE590ACA_Final_Project_Arthur\")

%Main Execution Script for Translational and Rotational Dynamics
%includes:
%  -Target and chaser state definitions
%  -Inertial Trajectory Propagation
%  -Relative Trajectory Propagation
%  -Various debugging
%  -Plotting
tic
clear all
close all

far_Rend = false;
unCoupled = false;  % Only Relevant when doing Close Range
% global u_history_global t_history_global
t_hist = [];


r_E = 6378*1000; %[m]
mu = 3.986004418e14;% [m3/s2]

%% propagate inertial target dynamics
% target
a_t = 11628*1000; %[m]
e_t = 0.4085;
%e_t = 0.2;%-----------debugging--------------
i_t = 70; % [deg]
%i_t = 89.5;%-----------debugging--------------
RAAN_t = 50; %[deg]
argp_t = 80; %[deg]
true_anom_t = 0; %[deg]

%set up sim time
t0 = 0;
tf = 8.5*3600; %[8.5 hours in sec] for Far Range (Default)
%alternatively plot for one orbit for now
% tf = 2*pi*sqrt(a_t^3/mu);
% tf = 50;
sim_tol = 1e-12;
dt = 30;
options = odeset('RelTol',sim_tol, 'AbsTol',sim_tol,'MaxStep',dt);
m = 1000;
B = [zeros(3,3); eye(3)]/m;

%define Q and R
Q = zeros(6,6);

Q(1:3,1:3) = eye(3).*0.01;
Q(4:6,4:6) = eye(3).*0.001;

R = eye(3).*10^4;

%original implementation: use the S_0 vec in the paper

r0 = [960; -590; 3290].*1000; %initial relative position, m
v0 = [0; -55; 0]; %initial relative velocity, m/s

DCM = rot_z(argp_t*180/pi)*rot_x(i_t*180/pi)*rot_x(RAAN_t*180/pi);
DCM = DCM';
rel_0 = [r0;v0];



q0          = [0;0;0;1];
if far_Rend
    % rel_0 = rel_0;
    h_wc0       = zeros(3,1);  % initial Reaction Wheel Output
    omega0      = [-0.3;0.5;0.1]*1e-2;  % Rad/sec
    % omega0      = [-0.3;0.5*0;0.1*0]*1e-2;  % Rad/sec
    % Assume that both body frames are aligned
    P1_c        = D(q0)*[1.5;1;0]*0;
    P0_t        = [1;0;1]*0;

else
    rel_0       = [25;10;50;0;-0.06;0];
    omega0      = [-0.4;0.5;0.2]*1e-2;  % Rad/sec
    % omega0      = -[-0.4;0.;0.]*1e-2;  % Rad/sec
    h_wc0       = [-3;5;1];  % For Near Term, use Figure 13!!
    h_wc0       = [ -1.49999846234356
        2.75000070984254
        0.600001631426673];  % From Far Range Sim
    % Assume that both body frames are aligned
    P1_c        = D(q0)*[1.5;1;0];
    P0_t        = [1;0;1];
    tf = 50;


    if unCoupled

        % If Pos/Att are uncoupled dynamics, we just set the relative state to be
        % the centered at the

        % Add Arbitrary docking points
        rel_0(1:3) =  rel_0(1:3) + P1_c - P0_t;  % Eqn (26)
        rel_0(4:6) =  rel_0(4:6) + cross(omega0,P1_c); % Eqn (27)

        P1_c        = D(q0)*[1.5;1;0]*0;
        P0_t        = [1;0;1]*0;
    end

    %% SCC Arbitrarily changing the Cost Weighting, this seems to yield the desired results
    POS_FACTOR = 100000;
    Q = POS_FACTOR*Q;
    R = R/POS_FACTOR;
end



tspan = linspace(t0,tf, 1000);

%target cartesian inertial position
[r_t0, v_t0] = keplerian2cartesian(a_t, e_t, i_t, RAAN_t, argp_t, true_anom_t, mu);
xt0 = [r_t0; v_t0];

[~, omegaT, ~] = kepler_orbital_elements_eval(0, mu, a_t, e_t);

omegaT_I = DCM*[0;0;omegaT];

xc0 = xt0 - [DCM*r0;DCM*(v0)-cross(omegaT_I ,r0)];  % Chaser

%propagate target with two body dynamics
inertial_t_trajectory = ode45(@(t,x) Cartesian_EOM(t,x,mu), tspan, xt0, options);
t_hist = deval(inertial_t_trajectory, tspan);
x_hist_t = t_hist(1:3,:);

%propagate chaser with two body dynamics
inertial_c_trajectory = ode45(@(t,x) Cartesian_EOM(t,x,mu), tspan, xc0, options);
c_hist = deval(inertial_c_trajectory, tspan);
x_hist_c = c_hist(1:3,:);

%


I_c = diag([500 550 600]);
I_t = I_c;
%
% t0 = 0;
% tf = 50;

x0 = [rel_0;q0;omega0;h_wc0];

% Add Arbitrary docking points
x0(1:3) =  x0(1:3) + P1_c - P0_t;  % Eqn (26)
x0(4:6) =  x0(4:6) + cross(omega0,P1_c); % Eqn (27)

%propagate with linearized relative orbital dynamics
%SDRE First
% relative_trajectory_SDRE = ode45(@(t,x) linearized_rel_orbital_dynamics(t, x, mu,a_t,e_t,'SDRE',B,R,Q), tspan, x0, options);
% state_hist_SDRE = deval(relative_trajectory_SDRE, tspan);

relative_trajectory_SDRE = ode45(@(t,x) coupledDynamics(t, x, mu,a_t,e_t,'SDRE',B,R,Q,I_c,P1_c,P0_t), tspan, x0, options);
state_hist_SDRE = deval(relative_trajectory_SDRE, tspan);

omega_hist = state_hist_SDRE(11:13,:)';
h_wc_hist = state_hist_SDRE(14:end,:)';
figure;
ylim_array = [-10 5;-5 10;-1 2]*1e-3;
axis_str = 'xyz';
for ii = 1:3
    subplot(3,1,ii)
    plot(tspan,omega_hist(:,ii));grid on
    ylabel(['$$\omega_{',axis_str(ii),'}$$[rad/sec]'],'Interpreter','latex')
    ylim(ylim_array(ii,:))
    set(gca,'TickLabelInterpreter','latex')
end
xlabel('time (sec)','Interpreter','latex')

figure;
ylim_array = [-4 0;0 6;0 2];
for ii = 1:3
    subplot(3,1,ii)
    plot(tspan,h_wc_hist(:,ii));grid on
    ylabel(['$$h_{WC_',axis_str(ii),'}$$ [N-m]'],'Interpreter','latex')
    % ylabel('N-m','Interpreter','latex')
    ylim(ylim_array(ii,:))
    set(gca,'TickLabelInterpreter','latex')
end
xlabel('time (sec)','Interpreter','latex')

%
% %LQR
% relative_trajectory_LQR = ode45(@(t,x) linearized_rel_orbital_dynamics(t, x, mu,a_t,e_t,'LQR',B,R,Q), tspan, x0, options);
% state_hist_LQR = deval(relative_trajectory_LQR, tspan);

%% plot relative state history
% figure()
% hold on
% plot3(state_hist_SDRE(1,:),state_hist_SDRE(2,:),state_hist_SDRE(3,:), 'r','LineWidth', 1.5)
% plot3(state_hist_SDRE(1,1),state_hist_SDRE(2,1),state_hist_SDRE(3,1), 'd','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
% plot3(state_hist_SDRE(1,end),state_hist_SDRE(2,end),state_hist_SDRE(3,end), 'd','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% legend('state history','initial', 'final')
% title('Relative State History')
% view(3)
if 1
    %relative position
    figure()
    sgtitle('Relative Position Time History')
    subplot(3,1,1)
    hold on
    plot(tspan, state_hist_SDRE(1,:), 'r','LineWidth', 1.5)
    %plot(tspan, state_hist_LQR(1,:), 'b--','LineWidth', 1.5)
    xlabel('time [sec]')
    ylabel('\Delta X (t)')
    legend('SDRE', 'LQR')

    subplot(3,1,2)
    hold on
    plot(tspan, state_hist_SDRE(2,:), 'r','LineWidth', 1.5)
    %plot(tspan, state_hist_LQR(2,:), 'b--','LineWidth', 1.5)
    xlabel('time [sec]')
    ylabel('\Delta Y (t)')
    legend('SDRE', 'LQR')

    subplot(3,1,3)
    hold on
    plot(tspan,state_hist_SDRE(3,:), 'r','LineWidth', 1.5)
    %plot(tspan, state_hist_LQR(3,:), 'b--','LineWidth', 1.5)
    xlabel('time [sec]')
    ylabel('\Delta Z (t)')
    legend('SDRE', 'LQR')


    %relative velocity
    figure()
    subplot(3,1,1)
    hold on
    sgtitle('Relative Velocity Time History')
    plot(tspan, state_hist_SDRE(4,:), 'r','LineWidth', 1.5)
    %plot(tspan, state_hist_LQR(4,:), 'b--','LineWidth', 1.5)
    xlabel('time [sec]')
    ylabel('\Delta v_x (t)')

    subplot(3,1,2)
    hold on
    plot(tspan, state_hist_SDRE(5,:), 'r','LineWidth', 1.5)
    %plot(tspan, state_hist_LQR(5,:), 'b--','LineWidth', 1.5)
    xlabel('time [sec]')
    ylabel('\Delta v_y (t)')

    subplot(3,1,3)
    hold on
    plot(tspan, state_hist_SDRE(6,:), 'r','LineWidth', 1.5)
    %  plot(tspan, state_hist_LQR(6,:), 'b--','LineWidth', 1.5)
    xlabel('time [sec]')
    ylabel('\Delta v_z (t)')


    %% Control history

    X_SDRE = state_hist_SDRE;
    % X_LQR = state_hist_LQR;

    for tt = 1:numel(tspan)

        t = tspan(tt);
        K_SDRE = calculateControllerGainfunction(t, X_SDRE(:,tt),  mu, a_t, e_t, 'SDRE',B,R,Q);
        % K_LQR = calculateControllerGainfunction(t, X_SDRE(:,tt),  mu, a_t, e_t, 'LQR',B,R,Q);
        u_history_global_SDRE(:,tt) = K_SDRE*X_SDRE(1:6,tt);
        %  u_history_global_LQR(:,tt) = K*X_LQR(:,tt);
    end

    t_history_global = tspan;

    %%
    figure()
    subplot(3,1,1)
    hold on
    sgtitle('Control History')
    plot(t_history_global, u_history_global_SDRE(1,:), 'r','LineWidth', 1.5)
    %   plot(t_history_global, u_history_global_LQR(1,:), 'b--','LineWidth', 1.5)
    xlabel('time [sec]')
    ylabel('u_x (t)')

    subplot(3,1,2)
    hold on
    plot(t_history_global, u_history_global_SDRE(2,:), 'r','LineWidth', 1.5)
    % plot(t_history_global, u_history_global_LQR(2,:), 'b--','LineWidth', 1.5)
    xlabel('time [sec]')
    ylabel('u_y (t)')

    subplot(3,1,3)
    hold on
    plot(t_history_global, u_history_global_SDRE(3,:), 'r','LineWidth', 1.5)
    % plot(t_history_global, u_history_global_LQR(3,:), 'b--','LineWidth', 1.5)
    xlabel('time [sec]')
    ylabel('u_z (t)')




end
toc



% Defining Primary Rotation Matrices
function DCM = rot_z(a)
% Author: Shaun Chemplavil
% Frame Rotation Matrix about 3rd Axis by angle 'a'
% input
% a: [deg]
% Output
% DCM: [3x3]  Direction Cosine Matrix, rotating frame about Z-axis

DCM = [cosd(a) sind(a) 0; -sind(a) cosd(a) 0; 0 0 1];
end

function DCM = rot_y(a)
% Author: Shaun Chemplavil
% Frame Rotation Matrix about 2nd Axis by angle 'a'
% input
% a: [deg]
% Output
% DCM: [3x3]  Direction Cosine Matrix, rotating frame about Y-axis
DCM = [cosd(a) 0 -sind(a); 0 1 0; sind(a) 0 cosd(a)];
end

function DCM = rot_x(a)
% Author: Shaun Chemplavil
% Frame Rotation Matrix about 1st Axis by angle 'a'
% input
% a: [deg]
% Output
% DCM: [3x3]  Direction Cosine Matrix, rotating frame about X-axis
DCM = [1 0 0;0 cosd(a) sind(a); 0 -sind(a) cosd(a)];
end
