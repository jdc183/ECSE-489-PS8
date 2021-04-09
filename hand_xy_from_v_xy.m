%function to integrate a velocity plan to derive the (x,y) plane
%input the velocity plan and the time step; returns the corresonding positions

function [hand_xy_plan] = hand_xy_from_v_xy (hand_xy_start,v_xy_plan,DT)
[dummy,npts]=size(v_xy_plan);
hand_xy_plan=zeros(2,npts);
hand_xy_plan(:,1)=hand_xy_start;
for i=2:npts
  hand_xy_plan(:,i)=hand_xy_plan(:,i-1)+0.5*(v_xy_plan(:,i-1)+v_xy_plan(:,i))*DT;
end

end
