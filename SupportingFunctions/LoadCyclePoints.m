function [RepClicks,ReplicationCluster]=LoadCyclePoints(M,ch,initval);
%Get cycle points and some first estimates from stored data inestead of
%clicking;
  RepClicks=M(ch).channels.RepClicks;
  [~,reps]=size(RepClicks);
  for rp=1:reps 
    rp
    ystart=RepClicks(rp).PosClick.firstframe;
    yend=RepClicks(rp).PosClick.lastframe;
    xstart=RepClicks(rp).PosClick.firstpos;
    xend=RepClicks(rp).PosClick.lastpos;
    poss=linspace(xstart,xend,yend-ystart+1);
    frs=[ystart:1:yend];
    ReplicationCluster(rp).PosKyTracCom.frames=frs;
    ReplicationCluster(rp).PosKyTracCom.clickpos=poss;
  end
end