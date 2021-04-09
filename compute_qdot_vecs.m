%given an array of hand velocity vectors and joint angles, compute the corresponding joint velocities

function [qdot_vecs] = compute_qdot_vecs (vxy_vecs,qvecs,a_vec)
  [dummy,npts]=size(vxy_vecs);
  qdot_vecs=zeros(2,npts);
  for isamp=1:npts
    J=compute_Jacobian(qvecs(:,isamp),a_vec);
    qdot_vec = inv(J)*vxy_vecs(:,isamp);
    qdot_vecs(:,isamp) = qdot_vec;
  end
end
