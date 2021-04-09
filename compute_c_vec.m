%given the h_ijk terms, and given qdot_vec, compute the centrifugal and coriolis torques
%SPECIALIZED FOR 2DOF PLANAR ARM
function [c_vec] = compute_c_vec (h_vals, qdot_vec)
    Cmat = zeros(2,2);
  Cmat(1,1) = (h_vals(1,2,2)+h_vals(1,2,1))*qdot_vec(2); 
  Cmat(1,2) = h_vals(1,2,2)*qdot_vec(2);
  Cmat(2,1) = h_vals(2,1,1)*qdot_vec(1);
  %Cmat(2,2)=0
  c_vec = Cmat*qdot_vec;

end
