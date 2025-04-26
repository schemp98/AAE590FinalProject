function dxdt = Cartesian_EOM(t,x,mu)

    r = x(1:3);
    v = x(4:6);
    norm_r = norm(r);
    
    f0 = [v; -mu/norm_r^3.*r];
   
    dxdt = f0;
return


