%computation of torques for 2DOF arm, given desired qddot_vec, qdot_vec and q_vec
%this routine incorporates the dynamic model constants

function [tau_vec,H,h_vals] = inv_dyn_2DOF (qddot_vec,qdot_vec,q_vec,DH_a_vec,grav_vec)
  %dynamic parameters:
  L1 = DH_a_vec(1); %1;
  Lc1=L1/2;
  L2 = DH_a_vec(2); %1;
  Lc2=L2/2;
  M1=1;
  M2=1;
  I1 = (1/12)*M1*L1*L1;
  I2 = (1/12)*M2*L2*L2;
  H = zeros(2,2);
  H(1,1) = M1*Lc1*Lc1+I1+M2*(L1*L1+Lc2*Lc2+2*L1*Lc2*cos(q_vec(2)))+I2;
  H(2,2) = M2*Lc2*Lc2+I2;
  H(1,2) = M2*L1*Lc2*cos(q_vec(2))+M2*Lc2*Lc2+I2;
  H(2,1) = H(1,2);
  %grav_vec=[0;0;9.8];
  h_vals=zeros(2,2,2);
  
  h111=0.0;
  h_vals(1,1,1) = h111;
  h122= -M2*L1*Lc2*sin(q_vec(2));
  h_vals(1,2,2)=h122;
  h112 = -2*M2*L1*Lc2*sin(q_vec(2));
  h_vals(1,1,2) = h112;
  h121=h112;
  h_vals(1,2,1)=h121;
  h211=M2*L1*Lc2*sin(q_vec(2));
  h_vals(2,1,1)=h211;
%  h222=0;
%  h212=0;
%  h221=0;
  
  
  [c_vec] = compute_c_vec (h_vals, qdot_vec);
  inertia_vec = H*qddot_vec;
  g_vec = grav_trqs(q_vec,grav_vec);
  tau_vec=inertia_vec+c_vec+g_vec;

end
