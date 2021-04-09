%given q_vec(t), compute max vx(t)

function [vx_max_vals] = compute_max_vels (q_vecs, qdot_max_vec,a_vec)
[dummy,npts]=size(q_vecs);
vx_max_vals=zeros(1,npts);
for i=1:npts
  J = compute_Jacobian(q_vecs(:,i),a_vec);
  constraint_vec = inv(J)*[1;0];
  vx_max_vals(i) = min(qdot_max_vec(1)/abs(constraint_vec(1)),qdot_max_vec(2)/abs(constraint_vec(2)));
  end
end
