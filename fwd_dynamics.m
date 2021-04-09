%fwd dynamics:

given q, qdot, (state) and  tau_vec (input, or excitation):
compute H(q)
compute C(q,qdot)
compute g(q)

tau_vec = H*qddot + C*qdot + g -->
qddot = H_inv*(tau_vec - C*qdot -g)
qdot(t+dt) = qdot(t)  + qddot*dt
q(t+dt) = q(t) + qdot(t)*dt


