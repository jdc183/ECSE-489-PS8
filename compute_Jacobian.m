%specialized for 2DOF planar arm:

function [J] = compute_Jacobian (q_vec, a_vec)
 %FK:  xy_vec(1) = a_vec(1)*cos(q_vec(1))+a_vec(2)*cos(q_vec(1)+q_vec(2));
 %     xy_vec(2) = a_vec(1)*sin(q_vec(1))+a_vec(2)*sin(q_vec(1)+q_vec(2));
 J=zeros(2,2);
 J(1,1) = -a_vec(1)*sin(q_vec(1))-a_vec(2)*sin(q_vec(1)+q_vec(2));  %dx/dq1
 J(1,2) = -a_vec(2)*sin(q_vec(1)+q_vec(2));%dx/dq2
 J(2,1) =   a_vec(1)*cos(q_vec(1))+a_vec(2)*cos(q_vec(1)+q_vec(2)); %dy/dq1
 J(2,2) =   a_vec(2)*cos(q_vec(1)+q_vec(2)); %dy/dq1

end
