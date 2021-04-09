%invent a proposed schedule of hand velocities; specify them at fixed intervals DT
%MUST start at 0 and end at 0

function [v_xy_plan] = plan_xy_velocities (DT,hand_start,hand_finish,q_dot_max,a_vec,tau_max_vec)
%choose npts:
hand_current = hand_start;
t = hand_finish - hand_start;
t = t/norm(t);
v_xy_plan = [];
n = 1;
vxy_prev = [0;0];%[1.311890684673591;0];
tau_vecs = [];
x_acc_min = [];
x_acc_max = [];
x_acc = [];
while norm(hand_current-hand_finish)>0.1
    
    q_vec = compute_IK(hand_current,a_vec);
    v_x_max = compute_max_vels(q_vec,q_dot_max,a_vec);
    vxy = [v_x_max;0];
    
    qdot_vec = compute_qdot_vecs(vxy,q_vec,a_vec);
    
    accxy = (vxy - vxy_prev)/DT;
    qddot_vec = compute_qddot_vecs(accxy,qdot_vec,q_vec,a_vec);
    [x_acc_min(n),x_acc_max(n)] = xddot_max(q_vec,qdot_vec,a_vec,tau_max_vec);
%     [tau,H,h_vals] = inv_dyn_2DOF(qddot_vec,qdot_vec,q_vec,a_vec,[0;0]);
%     if abs(tau(1)) > tau_max_vec(1) || abs(tau(2)) > tau_max_vec(2)
        if accxy(1) >  x_acc_max(n)
            accxy(1) = x_acc_max(n);
            vxy = vxy_prev + accxy*DT;
        elseif accxy(1) < x_acc_min(n)
            accxy(1) = x_acc_min(n);
            vxy = vxy_prev + accxy*DT;
        end
%     end
    hand_current = hand_current + vxy*DT;
    n = n+1
    x_acc(n) = accxy(1);
    qddot_vec = compute_qddot_vecs(accxy,q_vec,qdot_vec,a_vec,DT);
    [tau_vec,H,h_vals] = inv_dyn_2DOF (qddot_vec,qdot_vec,q_vec,a_vec,[0;-1]);
%     while (abs(tau_vec(1)) > abs(tau_max_vec(1)) || abs(tau_vec(2)) > abs(tau_max_vec(2)))
%         vxy = vxy_prev + accxy*DT*0.95;
%         qdot_vec = compute_qdot_vecs(vxy,q_vec,a_vec);
%         hand_current = hand_current + [v_x_max;0]*DT;
%         accxy = (vxy - vxy_prev)/DT;
%         qddot_vec = compute_qddot_vecs(accxy,q_vec,qdot_vec,a_vec,DT);
%         [tau_vec,H,h_vals] = inv_dyn_2DOF (qddot_vec,qdot_vec,q_vec,a_vec,[0;-1]);
%     end
    tau_vecs = [tau_vecs tau_vec];
    v_xy_plan = [v_xy_plan vxy];
    vxy_prev = vxy;
    
end


% %choose npts:
% hand_current = hand_finish;
% v_xy_plan = [];
% n = 1;
% vxy_prev = [0;0];%[1.311890684673591;0];
% tau_vecs = [];
% x_acc_min = [];
% x_acc_max = [];
% x_acc = [];
% % DT = DT;
% while norm(hand_current-hand_start)>0.1
%     
%     q_vec = compute_IK(hand_current,a_vec);
%     v_x_max = -compute_max_vels(q_vec,q_dot_max,a_vec);
%     vxy = [v_x_max;0];
%     
%     qdot_vec = compute_qdot_vecs(vxy,q_vec,a_vec);
%     
%     accxy = (vxy - vxy_prev)/DT;
%     qddot_vec = compute_qddot_vecs(accxy,qdot_vec,q_vec,a_vec);
%     [x_acc_min(n),x_acc_max(n)] = xddot_max(q_vec,qdot_vec,a_vec,tau_max_vec);
%     
%     if accxy(1) >  x_acc_max(n)
%         accxy(1) = x_acc_max(n);
%         vxy = vxy_prev + accxy*DT;
%     elseif accxy(1) < x_acc_min(n)
%         accxy(1) = x_acc_min(n);
%         vxy = vxy_prev + accxy*DT;
%     end
%     
%     hand_current = hand_current + vxy*DT;
%     n = n+1
%     x_acc(n) = accxy(1);
%     qddot_vec = compute_qddot_vecs(accxy,q_vec,qdot_vec,a_vec,DT);
%     [tau_vec,H,h_vals] = inv_dyn_2DOF (qddot_vec,qdot_vec,q_vec,a_vec,[0;-1]);
% 
%     tau_vecs = [tau_vecs tau_vec];
%     v_xy_plan = [v_xy_plan vxy];
%     vxy_prev = vxy;
%     
% end
% hand_xy_plan = hand_xy_from_v_xy (hand_start,v_xy_plan,DT);
figure(10)
plot(tau_vecs');
title('tau_vecs');
figure(11);
title('max/min accelerations')
plot(x_acc_max)
hold on
plot(x_acc_min)
plot(x_acc)
% figure
% vmax=0.2;
% dv = vmax/(npts/2);
% %create a plan for hand velocity as a function of time with timestep intervals DT
% v_xy_plan = zeros(2,npts);
% for i=2:npts/2
%   v_xy_plan(1,i) = v_xy_plan(1,i-1)+dv;
% end
% for i=1+npts/2:npts
%   v_xy_plan(1,i) = v_xy_plan(1,i-1)-dv;
% end
% %v_xy_plan(1,:)=1.0; %fixed velocity; not very clever
% v_xy_plan(1,1)=0.0; %required start/end velocity
% v_xy_plan(1,npts)=0.0;
end
