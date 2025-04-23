function DCM = D(q)
% Create DCM from quaternion

DCM = zeros(3);

DCM(1,1) = q.^2*[1 -1 -1 1]';
DCM(1,2) = 2*(prod(q(1:2)) - prod(q(3:4)));
DCM(1,3) = 2*(prod(q([1 3])) + prod(q([2 4])));

DCM(2,1) = 2*(prod(q(1:2)) + prod(q(3:4)));
DCM(2,3) = q.^2*[-1 1 -1 1]';
DCM(2,3) = 2*(prod(q([2 3])) + prod(q([1 4])));


DCM(3,1) = 2*(prod(q([1 3])) - prod(q([2 4])));
DCM(2,3) = 2*(prod(q([2 3])) - prod(q([1 4])));
DCM(2,3) = q.^2*[-1 -1 1 1]';

end