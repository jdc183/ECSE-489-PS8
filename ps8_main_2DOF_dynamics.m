%simulation of 2DOF arm dynamics
%coerce path along x=-1 to x=+1 at y=0.2
clear all
tau_vec_max = [10;5]
qdot_vec_max = [1;2]
a_vec = [1;1]; %link lengths
grav_vec=[0;0;9.8]; %gravity is parallel to joint axes, so no gravity torques
%specify the desired hand path: a line parallel to the x-axis, offset by y=0.2
x_start=-1;
x_end = 1;
y=0.2;
hand_xy_start = [x_start;y];
hand_xy_end = [x_end;y];
DT=0.001

%IMPORTANT: EDIT plan_xy_velocities() to put in a viable and desirable plan
%in terms of vx(t) (with vy(t)=0, for this example)
v_xy_plan = plan_xy_velocities(DT,hand_xy_start,hand_xy_end,qdot_vec_max,a_vec,tau_vec_max)
[dummy,npts]=size(v_xy_plan)
accel_xy_plan = accel_from_vel_plan(v_xy_plan,DT);
hand_xy_plan = hand_xy_from_v_xy(hand_xy_start,v_xy_plan,DT);
time_vec=0:DT:(npts-1)*DT;
figure(1)
plot(time_vec,hand_xy_plan(1,:),'r',time_vec,hand_xy_plan(2,:),'b');
title('planned x(t) (red) and y(t) (blue)')
xlabel('time (sec)')
ylabel('distance (m)')


figure(2)
plot(time_vec,v_xy_plan(1,:),'r',time_vec,v_xy_plan(2,:),'b');
title('planned vx(t) (red) and vy(t) (blue)')
xlabel('time (sec)')
ylabel('vel (m/s)')

figure(3)
plot(time_vec,accel_xy_plan(1,:),'r',time_vec,accel_xy_plan(2,:),'b');
title('planned accelerations: ax(t) (red) and ay(t) (blue)')
xlabel('time (sec)')
ylabel('accel (m/s/s)')



q_vecs=zeros(2,npts);
for isamp=1:npts
  hand_xy=hand_xy_plan(:,isamp);
  qvec = compute_IK(hand_xy,a_vec);
  q_vecs(:,isamp)=qvec;
end

figure(4)
plot(time_vec,q_vecs(1,:),'r',time_vec,q_vecs(2,:),'b')
title("IK: q1 (red), q2 (blue)")
xlabel('time (sec)')
ylabel('joint angle (rad)')

qdot_vecs = compute_qdot_vecs(v_xy_plan,q_vecs,a_vec);

figure(5)
clf
plot(time_vec,qdot_vecs(1,:),'r',time_vec,qdot_vecs(2,:),'b')
title('joint velocities: q1dot (red) and q2dot (blue)')
xlabel('time (sec)')
ylabel('joint velocity (rad/sec)')
hold on
plot(time_vec,qdot_vec_max(1)*ones(1,npts),'r')
plot(time_vec,-qdot_vec_max(1)*ones(1,npts),'r')
plot(time_vec,qdot_vec_max(2)*ones(1,npts),'b')
plot(time_vec,-qdot_vec_max(2)*ones(1,npts),'b')

%compute and plot corresponding joint accels:
qddot_vecs = compute_qddot_vecs (accel_xy_plan,q_vecs, qdot_vecs,a_vec,DT);


figure(6)
plot(time_vec,qddot_vecs(1,:),'r',time_vec,qddot_vecs(2,:),'b')
title('joint accelerations: q1ddot (red) and q2ddot (blue)')
xlabel('time (sec)')
ylabel('joint accel (rad/sec/sec)')

%compute the required torques for this plan:
tau_vecs = zeros(2,npts);
for i=1:npts
  %function [tau_vec,H,h_vals] = inv_dyn_2DOF (qddot_vec,qdot_vec,q_vec,DH_a_vec,grav_vec)
  [tau_vec,H,h_vals] = inv_dyn_2DOF(qddot_vecs(:,i),qdot_vecs(:,i),q_vecs(:,i),a_vec,grav_vec);
  tau_vecs(:,i)=tau_vec;
end

figure(7)
clf
plot(time_vec,tau_vecs(1,:),'r',time_vec,tau_vecs(2,:),'b')
title('joint torques: tau1 (red) and tau2 (blue)')
xlabel('time (sec)')
ylabel('joint toruqe (Nm/sec)')
hold on
plot(time_vec,tau_vec_max(1)*ones(1,npts),'r')
plot(time_vec,-tau_vec_max(1)*ones(1,npts),'r')
plot(time_vec,tau_vec_max(2)*ones(1,npts),'b')
plot(time_vec,-tau_vec_max(2)*ones(1,npts),'b')

[vx_max_vals] = compute_max_vels (q_vecs,qdot_vec_max,a_vec);
figure(8)
clf
plot(time_vec,vx_max_vals)
title('velocity constraints')
xlabel('time (sec)')
ylabel('vx max (m/sec)')

%test vel constraints:
q_test_vecs = zeros(2,npts);
for i=1:npts
  J=compute_Jacobian(q_vecs(:,i),a_vec);
  q_test_vecs(:,i) = inv(J)*[1;0]*vx_max_vals(i); 
end

figure(9)
plot(time_vec,q_test_vecs(1,:),time_vec,q_test_vecs(2,:))
title('test qdot max vals')
xlabel('time (sec)')
ylabel('qdot (rad/sec)')

return
% the rest of this is just debug/test code for supporting functions
%dx=0.01
%x_plan = x_start:dx:x_end;
%[dummy,npts]=size(x_plan)
%y_plan = y*ones(1,npts);
%hand_xy_samps = [x_plan;y_plan];

%TEST FK and IK:
%q_vec = [1;1]
%[xy_vec] = compute_FK (q_vec,a_vec)
%[q_vec] = compute_IK (xy_vec, a_vec)
%
%figure(1)
%plot(x_plan(1,:),y_plan(1,:))
%title('desired hand motion: x vs y')

%test the Jacobian:
dt=0.0001
dq1=0.001
dq2=0.001
for i=1:npts
  dq_test = [dq1;dq2];
  q_test = q_vecs(:,i)+dq_test; %qdot_vecs*dt;
  xy_nom = compute_FK(q_vecs(:,i),a_vec);
  xy_test = compute_FK (q_test,a_vec);
  dxy = xy_test - hand_xy_plan(:,i)
  J = compute_Jacobian (q_test, a_vec);
  Jxdq = J*dq_test
  %vxdt = vx_vec(i)*dt
end
%test use of J_inv:
%dx = [1;0]*vx*dt
for i=1:npts
  v_xy = v_xy_plan(:,i);
  q_nom = q_vecs(:,i); %qdot_vecs*dt;
  xy_nom = compute_FK(q_nom,a_vec);
  J = compute_Jacobian (q_test, a_vec);
  dq_vec = inverse(J)*v_xy*dt
  q_test = q_nom+dq_vec;
  q_err = dq_vec-qdot_vecs(:,i)*dt
  %pause
  xy_test = compute_FK (q_test,a_vec);
  dxy = xy_test - xy_nom

  %vxdt = vx_vec(i)*dt
end
%test compute Jacobian dot:
%qdot=[1;1];
%q=[1;1];
%dq=qdot*DT;
%J_nom = compute_Jacobian (q, a_vec);
%J_test = compute_Jacobian (q+dq, a_vec);
%J_dot = compute_Jacobian_dot (q, qdot,a_vec)  %compute_Jacobian_dot (q_vec, qdot_vec,a_vec)
%J_dot_approx = (1/DT)*(J_test-J_nom)




