
  function ThisDivInfo=Add_Forks(S,M,ch,ThisDivInfo)
  %Create artifical 'fork' traces spannning the diffusion phase; to be used
  %for setting domain borders for allocating division edges
      m=ThisDivInfo.family.me.idx;
      pp=ThisDivInfo.family.parent.idx;
      dL=ThisDivInfo.family.rightdaughter.idx;
      dR=ThisDivInfo.family.leftdaughter.idx;

      frs_m=S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames;
      xs_m= S(ch).channels.ReplicationCluster(m).PosKyTracCom.trackpos;
      plot(frs_m,xs_m, 'o-'); hold on;
      fr_MT=S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(end);   %mother termination frame
      x_MT=xs_m(end);  %mother termination pos
      if dL>0
         fr_LDI= S(ch).channels.ReplicationCluster(dL).PosKyTracCom.frames(1);     %Left daughter initiation frame 
         x_LDI=round(S(ch).channels.ReplicationCluster(dL).PosKyTracCom.trackpos(1));  % Left daughter initiation frame      
         frs_Lfork=fr_MT+1:fr_LDI-1;
         if length(frs_Lfork)==0,frs_Lfork=fr_MT;, end;      
         lf=length(frs_Lfork);     
         xs_Lfork=linspace(x_MT,x_LDI,lf);  
         Forks.left=[frs_Lfork' xs_Lfork' ];
         if lf>0
         plot(frs_Lfork,xs_Lfork,'ro-'); hold on;
         end
      else
          Forks.left=[fr_MT+1  x_MT];
          plot(fr_MT+1,  x_MT,'r*-'); hold on;
      end
     
      if dR>0
         fr_RDI= S(ch).channels.ReplicationCluster(dR).PosKyTracCom.frames(1);   %Right daughter initiation frame 
         x_RDI=round(S(ch).channels.ReplicationCluster(dR).PosKyTracCom.trackpos(1)); %Right daughter initiation frame
         frs_Rfork=fr_MT+1:fr_RDI-1;       
         if length(frs_Rfork)==0,frs_Rfork=fr_MT;, end;        
         lf=length(frs_Rfork);
         xs_Rfork=linspace(x_MT,x_RDI,lf);
         Forks.right=[frs_Rfork' xs_Rfork' ];  
         if lf>0
         plot(frs_Rfork,xs_Rfork,'mo-'); hold on;
         end
      else
          Forks.right=[fr_MT+1  x_MT];
          plot(fr_MT+1,  x_MT,'mo-'); hold on;
      end
      ThisDivInfo.Forks=Forks;
  end   