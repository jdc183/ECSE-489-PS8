%given a desired time-history of x and y velocities and a uniform time step, deduce the corresponding x and y accelerations
function [accel_xy_plan] = accel_from_vel_plan (v_xy_plan,dt)
  [dummy,npts] = size(v_xy_plan);
  accel_xy_plan=zeros(2,npts);
  for i=1:npts-1
    accel_xy_plan(:,i)=(v_xy_plan(:,i+1)-v_xy_plan(:,i))/dt;
  end
  

end
