function B = computeControl(q,I_C)

% eqn (50)
GO = -D(q)*inv(I_C);  %  I'm not sure what the sign should be.....

B = zeros(6,3);
B(4:6,:) = GO;