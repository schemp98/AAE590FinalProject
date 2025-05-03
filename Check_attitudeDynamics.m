
close all
clear all
%% Test out attitude dynamics

att0 = zeros(3,1);

q0 = [0;0;0;1];
omega0 = [-0.4;0.5;0.2];  % Rad/sec
% omega0 = [-0.0;0.5;0.];  % Rad/sec

I_c = diag([500 550 600]);
I_t = I_c;

t0 = 0;
tf = 50;

x0 = [q0;omega0];
u = [0.4*3;0;0.6];


dt = 1e-3;

tol = 1e-12;
odeOptions = odeset('RelTol',tol,'AbsTol',tol);

out = ode45(@attitudeDynamics,[t0 tf],x0,odeOptions,I_c,u);

plot(out.x,out.y)

P = [1;0;0];

q = out.y(1:4,:);

%%
for ii = 1:numel(out.x)

    qq = q(:,ii);
    P_array(:,ii) = D(qq)*P;

end


figure;plot(out.x,P_array');legend('x','y','z')
max_check = abs(max(vecnorm(P_array)) - norm(P) ) < 1e-10;  % Check that our vector did not change size
min_check = abs(min(vecnorm(P_array)) - norm(P) ) < 1e-10;  % Check that our vector did not change size

if not(or(max_check,min_check))
error('something went wrong')
end