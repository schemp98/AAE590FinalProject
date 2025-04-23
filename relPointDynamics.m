function dx = relPointDynamics(t,x,P_C,P_T,w)

% Eqn (31)
P1 = [w(2)*(w(1)*P_C(2) - w(2)*P_C(1))-w(3)*(w(3)*P_C(1)-w(1)*P_C(3))]...
    +wd(2)*P_C(3)-wd(3)*P_C(2) + 2*wOT*(-w(3)*P_C(1) + w(1)*P_C(3))...
    +wdOT*(-P_C(2)+P_T(2))+3*wOT^2*(-P_C(1)+P_T(1)) + mu/r_t^2 - mu*r_t/r_c^3;

P2 = [w(3)*(w(2)*P_C(3) - w(3)*P_C(2))-w(1)*(w(1)*P_C(2)-w(2)*P_C(1))]...
    +wd(3)*P_C(1)-wd(1)*P_C(3) + 2*wOT*(-w(2)*P_C(3) + w(3)*P_C(2))...
    +wdOT*(-P_C(1)+P_T(1))+wOT^2*(-P_C(2)+P_T(2));


P3 = [w(1)*(w(3)*P_C(1) - w(1)*P_C(3))-w(2)*(w(2)*P_C(3)-w(3)*P_C(2))]...
    +wd(1)*P_C(2)-wd(2)*P_C(1) + wOT*(-P_C(3)+P_T(3));



dx(1:3) = x(4:6);

dx(4)   =  2*wOT*x(5) + wdOT*x(2) + 3*wOT^2*x(1);
dx(5)   = -2*wOT*x(4) - wdOT*x(1);
dx(6)   = -wOT^2*x(3) - wdOT*x(1);

 dx(4:6) = dx(4:6) + a + [P1;P2;P3];