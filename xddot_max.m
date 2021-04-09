%compute max accel along a path subject to torque saturation
%this function is not completed--only outlined

function [x_acc_min,x_acc_max] = xddot_max (q_vec, qdot_vec,a_vec,tau_max_vec)
% global Q QDOT A TAU_MAX
% Q = q_vec;
% QDOT = qdot_vec;
% A = a_vec;
% TAU_MAX = tau_max_vec;
% 
% x_acc_max = abs(fminunc(@tau_max_err,1));
% x_acc_min = -x_acc_max;
%tau = H qddot + C qdot
%qddot = J_inv*vdot - J_inv*Jdot*qdot
%tau = -H*J_inv*Jdot*qdot + C*qdot + H*J_inv*[1;0]*ax
% axy = (H*J_inv)\(tau + H*J_inv*Jdot*qdot - C*qdot)
% ax = ax(1)
%find min/max ax subject to tau_max constraints

g_vec = [0;0];%grav_trqs (q_vec,[0;-1]);
[~,H,h_vals] = inv_dyn_2DOF ([1;1],qdot_vec,q_vec,a_vec,g_vec);
[c_vec] = compute_c_vec (h_vals, qdot_vec);

J = compute_Jacobian(q_vec, a_vec);
Jdot = compute_Jacobian_dot(q_vec,qdot_vec, a_vec);

upper = 1;
lower = -1;
% tau = -H*J\Jdot*qdot_vec + c_vec + H*J\[1;0]*upper;
[tau,~,~] = inv_dyn_2DOF (compute_qddot_vecs ([upper;0],q_vec,qdot_vec,a_vec,.001),qdot_vec,q_vec,a_vec,g_vec);

while abs(tau(1))<abs(tau_max_vec(1)) || abs(tau(2))<abs(tau_max_vec(2))
     lower = upper;
    upper = 2*upper;
%     tau = -H*J\Jdot*qdot_vec + c_vec + H*J\[1;0]*upper;
    qddot_vec = compute_qddot_vecs ([upper;0],q_vec,qdot_vec,a_vec,.001);
    [tau,~,~] = inv_dyn_2DOF (qddot_vec,qdot_vec,q_vec,a_vec,g_vec);
end
% while abs(tau(1))<abs(tau_max_vec(1)) || abs(tau(2))<abs(tau_max_vec(2))
%     upper = lower; 
%     lower = 2*lower;
% %     tau = -H*J\Jdot*qdot_vec + c_vec + H*J\[1;0]*upper;
%     [tau,~,~] = inv_dyn_2DOF (compute_qddot_vecs ([lower;0],q_vec,qdot_vec,a_vec,.001),qdot_vec,q_vec,a_vec,g_vec);
% end
err = 8;
i = 0;
test = 0;
while err > 0.01 && i < 50
    test = (upper + lower)/2;
%     upper = upper;
%     lower = lower;%J\([test;0]-Jdot*qdot_vec)
%     tau = -H*J\Jdot*qdot_vec + c_vec + H*J\[1;0]*test
    qddot_vec = compute_qddot_vecs ([upper;0],q_vec,qdot_vec,a_vec,.001);
    [tau,~,~] = inv_dyn_2DOF (qddot_vec,qdot_vec,q_vec,a_vec,g_vec);
    err = min(abs(abs(tau) - abs(tau_max_vec)));
    i = i+1;
    if abs(tau(1))>abs(tau_max_vec(1)) || abs(tau(2))>abs(tau_max_vec(2))
        upper = test;
    else
        lower = test;
    end
end

x_acc_max = min(11,test*.95);
x_acc_min = max(-11,-test*.95);

end

% function err = tau_max_err(x_acc)
%     global Q QDOT A TAU_MAX
%     qddot_vec = compute_qddot_vecs ([x_acc;0],Q,QDOT,A,.001);
%     [tau,~,~] = inv_dyn_2DOF (qddot_vec,QDOT,Q,A,[0;0]);
%     err = min(abs(abs(tau) - abs(TAU_MAX)));
% end
% % tau_max_vec = dot(tau,tau_max_vec)*tau;
% 
% 
% qddot_vec_1 = H\(tau_max_vec - c_vec - g_vec);
% a_max_1 = J*qddot_vec_1 + Jdot*qdot_vec;
% 
% qddot_con = 
% % % a_max_1 = (H/J)\(tau_max_vec + H*J\Jdot*qdot_vec - c_vec - g_vec);
% % 
% % qddot_vec_2 = H\(-tau_max_vec - c_vec - g_vec);
% % a_max_2 = J*qddot_vec_2 + Jdot*qdot_vec;
% % % a_max_2 = (H/J)\(-tau_max_vec + H*J\Jdot*qdot_vec - c_vec - g_vec);
% % 
% % qddot_vec_3 = H\([1 0;0 -1]*tau_max_vec - c_vec - g_vec);
% % a_max_3 = J*qddot_vec_3 + Jdot*qdot_vec;
% % % a_max_3 = (H/J)\([1 0;0 -1]*tau_max_vec + H*J\Jdot*qdot_vec - c_vec - g_vec);
% % 
% % qddot_vec_4 = H\(-[1 0;0 -1]*tau_max_vec - c_vec - g_vec);
% % a_max_4 = J*qddot_vec_4 + Jdot*qdot_vec;
% % % a_max_4 = (H/J)\(-[1 0;0 -1]*tau_max_vec + H*J\Jdot*qdot_vec - c_vec - g_vec);
% % 
% % x_acc_max = max([a_max_1(1) a_max_2(1) a_max_3(1) a_max_4(1)]);
% % x_acc_min = min([a_max_1(1) a_max_2(1) a_max_3(1) a_max_4(1)]);
% x_acc_max = lower;
% x_acc_min = -lower;
% end
