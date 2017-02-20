function A060_RepliCluster00_Click(user,exp,WorkspaceOutname)

  if nargin<2
    exp='Exp001_DnaN_TUS_dif_01092016_M';
  end
  if nargin<1
    user = 'MarkPC';
  end

  %Analyze_ReplicationCluster
  %ReplicationCluster
  %-------------------------------------------------------------------------
  tic
  close all;

  % exp='001_DnaN_TUS_dif_30122014_difsignal';

  actions.init=1;        %default=1 (new analysis)
  actions.clickcycles=1; %default=1 (new analysis) else load
  actions.divisionoverlay=0;
  actions.eliminate_growth=1;

  peak_tol=0;   %sensitivity for detecting peaks

  initval=A001_Images_Set_Experiment(user,exp);

  if nargin == 3;
    initval.nms = WorkspaceOutname;
  end

  chans=initval.channelno;
  for ch=1:chans
    close all

    fprintf('Channel %d of %d.\n', ch, chans);

    DnaNIdx=find(ismember(initval.viewchan,initval.DnaNchan));
    Channelpath=char(strcat(initval.basepath,initval.nms{ch}{DnaNIdx},'.mat'));
    load(Channelpath, 'chanstk_BF','chanstk_FL','endpoints', 'kymo_BF','kymo_FL','presets');
    [r,c]=size(kymo_BF); 

    if actions.divisionoverlay
      BW=MakeBinaryEdgeImagefrom(kymo_BF,peak_tol);
      kymo_FL_click(:,:,1)=kymo_FL*255/range(kymo_FL(:))+0.5*BW*255;
    else
      kymo_FL_click=kymo_FL;
    end

    %static_kymo=Processing_Straighten_Growth(kymo_FL,initval);  %TEST

    kymoprops.width=c;
    kymoprops.duration=r;
    kymoprops.zoom=initval.zoom;  %used for clicking

    if actions.init
      [ReplicationCluster,RepClicks]=RepliCluster_Init(ch,kymo_FL_click,kymoprops, initval,actions);  %Init start spots
    end
    if actions.clickcycles
      figure(1)
      [ReplicationCluster,RepClicks]=RepliCluster_GetCyclePoints(ch,ReplicationCluster,RepClicks,kymoprops, kymo_FL_click,initval,actions);
      %Detect (manually) start and end points of a replication spot 'cluster';
      %that may consist of one or two spots but is considered to represent one
      %replication cycle

    else
      %load  database--------------------------------------------------
      outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
      load(outname_usr,'M');
      [RepClicks,ReplicationCluster]=LoadCyclePoints(M,ch,initval);
    end

    ChanNum=size(initval.viewchan,2);

    for i=1:ChanNum
      kymoprops.WorkspaceOutName=char(initval.nms{ch}{i}); 
      outname=strcat(initval.basepath,kymoprops.WorkspaceOutName);
      save(outname, 'initval', 'RepClicks', 'ReplicationCluster',  '-append');
    end
  end 
  fprintf('A060 has completed.\n');
end