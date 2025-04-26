function E = mean_to_ecc_anomaly_2(M, e)
tol = 1e-9;
E = M; % initial guess
while abs(M - (E - e*sin(E))) > tol
    E = E - (E - e*sin(E) - M) / (1 - e*cos(E));
end

E = mod(E, 2*pi);
end

