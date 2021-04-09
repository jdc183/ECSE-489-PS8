%specialized for 2DOF planar arm: given desired (x,y) and link lengths in a_vec, compute
%the "elbow up" solution for IK and return this in q_vec
%uses soln from text, 3-DOF planar arm IK, pg 68

function [q_vec] = compute_IK (hand_xy, a_vec)
  q_vec=zeros(2,1);
  c2 = (hand_xy(1)*hand_xy(1)+hand_xy(2)*hand_xy(2)-a_vec(1)*a_vec(1)-a_vec(2)*a_vec(2))/(2*a_vec(1)*a_vec(2));
  q_vec(2) = acos(c2);
  s2 = sin(q_vec(2));
  r_sqd = hand_xy'*hand_xy;
  s1 = ((a_vec(1)+a_vec(2)*c2)*hand_xy(2)-a_vec(2)*s2*hand_xy(1))/r_sqd;
  c1 = ((a_vec(1)+a_vec(2)*c2)*hand_xy(1)+a_vec(2)*s2*hand_xy(2))/r_sqd;
  q_vec(1) = atan2(s1,c1);

end
