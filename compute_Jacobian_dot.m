%compute_Jacobian_dot: given q_vec and qdot_vec, compute d/dt(Jacobian)

function [Jdot] = compute_Jacobian_dot (q_vec, qdot_vec,a_vec)
 Jdot=zeros(2,2);
 %J(1,1) = -a_vec(1)*sin(q_vec(1))-a_vec(2)*sin(q_vec(1)+q_vec(2));  %dx/dq1
 %J(1,2) = -a_vec(2)*sin(q_vec(1)+q_vec(2));%dx/dq2
 %J(2,1) =   a_vec(1)*cos(q_vec(1))+a_vec(2)*cos(q_vec(1)+q_vec(2)); %dy/dq1
 %J(2,2) =   a_vec(2)*cos(q_vec(1)+q_vec(2)); %dy/dq1
 Jdot(1,1) = -a_vec(1)*cos(q_vec(1))*qdot_vec(1)-a_vec(2)*cos(q_vec(1)+q_vec(2))*(qdot_vec(1)+qdot_vec(2));
 Jdot(1,2) = -a_vec(2)*cos(q_vec(1)+q_vec(2))*(qdot_vec(1)+qdot_vec(2));
 Jdot(2,1) =  -a_vec(1)*sin(q_vec(1))*qdot_vec(1)-a_vec(2)*sin(q_vec(1)+q_vec(2))*(qdot_vec(1)+qdot_vec(2));
 Jdot(2,2) = -a_vec(2)*sin(q_vec(1)+q_vec(2))*(qdot_vec(1)+qdot_vec(2));
end
