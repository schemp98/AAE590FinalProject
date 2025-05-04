function DCM = D(q)
% Create DCM from quaternion Eqn (11)

DCM = zeros(3);

DCM(1,1) = [1 -1 -1 1]*q.^2;
DCM(1,2) = 2*(prod(q(1:2)) - prod(q(3:4)));
DCM(1,3) = 2*(prod(q([1 3])) + prod(q([2 4])));

DCM(2,1) = 2*(prod(q(1:2)) + prod(q(3:4)));
DCM(2,2) = [-1 1 -1 1]*q.^2;
DCM(2,3) = 2*(prod(q([2 3])) - prod(q([1 4])));


DCM(3,1) = 2*(prod(q([1 3])) - prod(q([2 4])));
DCM(3,2) = 2*(prod(q([2 3])) + prod(q([1 4])));
DCM(3,3) = [-1 -1 1 1]*q.^2;

DCM = DCM';  % Equation is incorrect... this should be transposed
end