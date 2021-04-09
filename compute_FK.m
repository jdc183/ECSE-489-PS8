%compute FK, specialized for planar 2-link arm

function [xy_vec] = compute_FK (q_vec,a_vec)
 xy_vec=zeros(2,1);
 xy_vec(1) = a_vec(1)*cos(q_vec(1))+a_vec(2)*cos(q_vec(1)+q_vec(2));
 xy_vec(2) = a_vec(1)*sin(q_vec(1))+a_vec(2)*sin(q_vec(1)+q_vec(2));
end
