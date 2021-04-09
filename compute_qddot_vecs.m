%compute joint accelerations, given q_vecs and qdot_vecs
%compute the joint accelerations, qddot_vecs
%follows from twist = J*qdot
%d/dt(twist) = J*qddot + Jdot*qdot
%--> J_inv*( d/dt(twist) - Jdot*qdot) = qddot
function [qddot_vecs] = compute_qddot_vecs (accel_xy_vecs,q_vecs, qdot_vecs,a_vec,DT)
[dummy,npts]=size(q_vecs);
qddot_vecs=zeros(2,npts);
for i=1:npts
  J = compute_Jacobian(q_vecs(:,i),a_vec);
  Jdot = compute_Jacobian_dot(q_vecs(:,i),qdot_vecs(:,i),a_vec);
  qddot_vecs(:,i) = J\(accel_xy_vecs(:,i)-Jdot*qdot_vecs(:,i));
end

end
