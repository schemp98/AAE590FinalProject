%Main Execution Scrip for Translational Dynamics
%includes:
%  -Target and chaser state definitions
%  -Inertial Trajectory Propagation
%  -Relative Trajectory Propagation
%  -Various debugging
%  -Plotting
tic
clear all
close all

% global u_history_global t_history_global
t_hist = [];
% global drdt_log t_log omegaT_log omegaT_dot_log
% drdt_log = [];



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
%tf = 8.5*3600; %[8.5 hours in sec]
%alternatively plot for one orbit
tf = 2*pi*sqrt(a_t^3/mu);
sim_tol = 1e-12;
dt = 30;
options = odeset('RelTol',sim_tol, 'AbsTol',sim_tol,'MaxStep',dt);
%options = odeset('RelTol',sim_tol, 'AbsTol',sim_tol);

B = [zeros(3,3); eye(3)];

%define Q and R
Q = zeros(6,6);

Q(1:3,1:3) = eye(3).*0.01;
Q(4:6,4:6) = eye(3).*0.001;

R = eye(3).*10^4;


tspan = linspace(t0,tf, 3000);

%target cartesian inertial position
[r_t0, v_t0] = keplerian2cartesian(a_t, e_t, i_t, RAAN_t, argp_t, true_anom_t, mu)
xt0 = [r_t0; v_t0];


%instead, define chaser through keplerian elements to debug:
% chaser
a_c = 11628*1000; %[m]
e_c = 0.4085;
%e_c = 0.6;%-----------debugging--------------
i_c = 70 + 5; % [deg]
%i_c = 89.5;%-----------debugging--------------
RAAN_c = 50; %[deg]
argp_c = 80; %[deg]
true_anom_c = 0.1; %[deg]

%target cartesian inertial position
[r_c0, v_c0] = keplerian2cartesian(a_c, e_c, i_c, RAAN_c, argp_c, true_anom_c, mu);
xc0 = [r_c0; v_c0];

rel_0 = xt0 -xc0;

%original implementation: use the S_0 vec in the paper
r0 = [960; -590; 290].*1000; %initial relative position, m
v0 = [0; -55; 0]; %initial relative velocity, m/s
rel_0 = [r0;v0];
xc0 = xt0 - rel_0;



%propagate target with two body dynamics
inertial_t_trajectory = ode45(@(t,x) Cartesian_EOM(t,x,mu), tspan, xt0, options);
t_hist = deval(inertial_t_trajectory, tspan);
x_hist_t = t_hist(1:3,:);

%propagate chaser with two body dynamics
inertial_c_trajectory = ode45(@(t,x) Cartesian_EOM(t,x,mu), tspan, xc0, options);
c_hist = deval(inertial_c_trajectory, tspan);
x_hist_c = c_hist(1:3,:);


%plot target orbit
figure()
hold on
plot3(x_hist_t(1,:),x_hist_t(2,:),x_hist_t(3,:), 'r','LineWidth', 1.5)
plot3(x_hist_c(1,:),x_hist_c(2,:),x_hist_c(3,:), 'b--','LineWidth', 1.5)
plot3(x_hist_c(1,1),x_hist_c(2,1),x_hist_c(3,1), 'd','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10)
plot3(x_hist_t(1,1),x_hist_t(2,1),x_hist_t(3,1), 'd','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
[earthX, earthY, earthZ] = sphere;
surf(r_E*earthX,r_E*earthY,r_E*earthZ, 'FaceColor','none')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
title('Inertially Propagated Orbits')
legend('target', 'chaser')
hold off

%plot difference in inertial frame
figure()
hold on
plot3(x_hist_t(1,:)-x_hist_c(1,:),x_hist_t(2,:)-x_hist_c(2,:),x_hist_t(3,:)-x_hist_c(3,:), 'r','LineWidth', 1.5)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('Difference of Inertially Propagated Orbits')
view(3)
hold off
%
% for i = 1:size(x_hist_t,2)
%     normdiff(i) = norm(x_hist_t(:,i)-x_hist_c(:,i));
% end
% %plot magnitude of difference in inertial frame
% figure()
% hold on
% plot(normdiff-  normdiff(1))
% xlabel('t')
% ylabel('norm(diff)')
% title('Magnitude of Difference of Inertially Propagated Orbits')
% hold off
%% propagate relative dynamics

x0 = rel_0;

%propagate with linearized relative orbital dynamics
%SDRE First
relative_trajectory_SDRE = ode45(@(t,x) linearized_rel_orbital_dynamics(t, x, mu,a_t,e_t,'SDRE',B,R,Q), tspan, x0, options);
state_hist_SDRE = deval(relative_trajectory_SDRE, tspan);

%LQR
relative_trajectory_LQR = ode45(@(t,x) linearized_rel_orbital_dynamics(t, x, mu,a_t,e_t,'LQR',B,R,Q), tspan, x0, options);
state_hist_LQR = deval(relative_trajectory_LQR, tspan);

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

%relative position
figure()
sgtitle('Relative Position Time History')
subplot(3,1,1)
hold on
plot(tspan, state_hist_SDRE(1,:), 'r','LineWidth', 1.5)
plot(tspan, state_hist_LQR(1,:), 'b--','LineWidth', 1.5)
xlabel('time [sec]')
ylabel('\Delta X (t)')
xlim([0 500])
legend('SDRE', 'LQR')
subplot(3,1,2)
hold on
plot(tspan, state_hist_SDRE(2,:), 'r','LineWidth', 1.5)
plot(tspan, state_hist_LQR(2,:), 'b--','LineWidth', 1.5)
xlabel('time [sec]')
ylabel('\Delta Y (t)')
legend('SDRE', 'LQR')
xlim([0 500])
subplot(3,1,3)
hold on
plot(tspan,state_hist_SDRE(3,:), 'r','LineWidth', 1.5)
plot(tspan, state_hist_LQR(3,:), 'b--','LineWidth', 1.5)
xlabel('time [sec]')
ylabel('\Delta Z (t)')
legend('SDRE', 'LQR')
xlim([0 500])

%relative velocity
figure()
subplot(3,1,1)
hold on
sgtitle('Relative Velocity Time History')
plot(tspan, state_hist_SDRE(4,:), 'r','LineWidth', 1.5)
plot(tspan, state_hist_LQR(4,:), 'b--','LineWidth', 1.5)
xlabel('time [sec]')
ylabel('\Delta v_x (t)')
xlim([0 500])
subplot(3,1,2)
hold on
plot(tspan, state_hist_SDRE(5,:), 'r','LineWidth', 1.5)
plot(tspan, state_hist_LQR(5,:), 'b--','LineWidth', 1.5)
xlabel('time [sec]')
ylabel('\Delta v_y (t)')
xlim([0 500])
subplot(3,1,3)
hold on
plot(tspan, state_hist_SDRE(6,:), 'r','LineWidth', 1.5)
plot(tspan, state_hist_LQR(6,:), 'b--','LineWidth', 1.5)
xlabel('time [sec]')
ylabel('\Delta v_z (t)')
xlim([0 500])

%% Control history

X = state_hist_SDRE;

for tt = 1:numel(tspan)

t = tspan(tt);
K = calculateControllerGainfunction(t,  mu, a_t, e_t, 'SDRE',B,R,Q);
u_history_global(:,tt) = K*X(:,tt);
end

t_history_global = tspan;

%%
figure()
subplot(3,1,1)
sgtitle('Control History')
plot(t_history_global, u_history_global(1,:), 'g','LineWidth', 1.5)
xlabel('time [sec]')
ylabel('u_x (t)')
xlim([0 500])
subplot(3,1,2)
plot(t_history_global, u_history_global(2,:), 'g','LineWidth', 1.5)
xlabel('time [sec]')
ylabel('u_y (t)')
xlim([0 500])
subplot(3,1,3)
plot(t_history_global, u_history_global(3,:), 'g','LineWidth', 1.5)
xlabel('time [sec]')
ylabel('u_z (t)')
xlim([0 500])



%% debugging and other stuff, not used for replicating results

% %% checking dr/dt time history
% figure()
% plot(t_log,drdt_log)
% title('drdt time history')
%
% figure()
% plot(t_log,omegaT_dot_log)
% title('\dot{\omega_T} time history')
%
% omegaT_dot_grad= gradient(omegaT_log', t_log)';
% figure()
% plot(t_log,omegaT_dot_grad)
% title('\omega_T gradient time history')
%
% figure()
% plot(t_log,(omegaT_dot_log -omegaT_dot_grad))
% title('\omega_T gradient vs \dot{\omega_T} difference time history')

% %% testing - reconstruct chaser trajectory in inertial frame
% rT_hist = t_hist(1:3,:)';
% vT_hist = t_hist(4:6,:)';
%
% rC_hist = c_hist(1:3,:)';
% vC_hist = c_hist(4:6,:)';
%
% x_rel_hist = state_hist';
%
% [rC_hist, vC_hist] = reconstruct_chaser_from_relative(tspan, x_rel_hist, rT_hist, vT_hist);
%
% figure()
% hold on
% plot3(rC_hist(:,1),rC_hist(:,2),rC_hist(:,3), 'r','LineWidth', 1.5)
% plot3(rT_hist(:,1),rT_hist(:,2),rT_hist(:,3), 'b--','LineWidth', 1.5)
% plot3(rC_hist(1,1),rC_hist(1,2),rC_hist(1,3), 'd','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
% plot3(rT_hist(1,1),rT_hist(1,2),rT_hist(1,3), 'd','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10)
% [earthX, earthY, earthZ] = sphere;
% surf(r_E*earthX,r_E*earthY,r_E*earthZ, 'FaceColor','none')
% view(3)
% axis equal
% xlabel('x')
% ylabel('y')
% zlabel('z')
% xlim([-2*r_E, 2*r_E])
% ylim([-2*r_E, 2*r_E])
% zlim([-2*r_E, 2*r_E])
%plot3(chaser_state(1,1),chaser_state(2,1),chaser_state(3,1), 'd','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)

% %%
% rel_LVLH_hist = inertial_diff_to_LVLH(rT_hist, vT_hist, rC_hist);
%
% figure;
% plot3(rel_LVLH_hist(:,1), rel_LVLH_hist(:,2), rel_LVLH_hist(:,3));
% xlabel('x_{LVLH} [m]');
% ylabel('y_{LVLH} [m]');
% zlabel('z_{LVLH} [m]');
% title('Relative Motion in Target LVLH Frame');
% grid on;
%
% view(3);

%reconstruct control history:
%using global var for now

toc

