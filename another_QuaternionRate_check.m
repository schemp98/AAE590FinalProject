%% DCM Check
q0          = [0;0;0;1];

w0 = zeros(3,1);
% w0(randi([1 3],1)) = randn*1e-3
w0(1) = randn*1e-3

% w0 = randn(3,1)
t = 1;
qd = relQuatRate(t,q0,w0);

q = q0 + qd; q = q/norm(q);
DCM = D(q)
% Use small angle approximation
phi   = DCM(2,3);  % roll
theta = DCM(3,1);  % pitch
psi   = DCM(1,2);  % yaw




% Eqn (12)
QQ =@(q) [ q(4)*eye(3) + [0 -q(3) q(2);q(3) 0 -q(1);-q(2) q(1) 0];-q(1:3)'];

qd - 0.5*QQ(q0)*w0

PP = [1;1;0]

PP_paper = DCM*PP
DCM_scc = rot_x(phi*180/pi)*rot_y(theta*180/pi)*rot_z(psi*180/pi);
PP_scc = DCM_scc*PP



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
