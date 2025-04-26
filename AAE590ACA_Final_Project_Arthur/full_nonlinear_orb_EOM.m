% function dx = full_nonlinear_orb_EOM(t,x,mu,r_c,u)
% 
% 
%     %implement dynamics per slide 41 (no need for rcdot, circular orbit)
%     rd = [r_c + x(1); x(2); x(3)];
%     norm_rd = norm(rd);
% 
%     n = sqrt(mu/r_c^3);
% 
%     A = [zeros([3 3]) eye(3); n^2-mu/norm_rd^3 0 0 0 2*n 0; ...
%         0 n^2-mu/norm_rd^3 0 -2*n 0 0; 0 0 -mu/norm_rd^3 0 0 0];
% 
%     B = [zeros([3 3]); eye(3)];
% 
%     c = [0; 0; 0; mu/r_c^2-mu/norm_rd^3*r_c; 0; 0];
% 
% 
% 
%     %fix the last index - one more state than control input
%     if idx < 51
%         dx = A*x + B*u(:,idx) +c;
%     else
%          dx = A*x + B*[0;0;0] +c;
%     end
% 
% 
% 
% return