%function A102_Processing_Find_Division_Times(exp)
function Processing_Find_Division_Times(exp,user,ColourIdx)
  %Load database, 
  %1) For every replicationcluster....
  %find the two daugthers
  %(if found) create 'left'and 'right' forks that start at termination and extend to the new
  %daugther init points. The new division should lie in the region encompassed by this fork and the subsequent region.
  %this division marks the end of this bacterium.
  %We use a binary brightfield-kymograph to do this

  %JacobKers 2012

  close all
  actions.savedata=1;         %default=1 (new analysis)
  actions.loaddatabase=1;     %default=1 (new analysis)
  actions.plot=1;             %default=0 (quick run)

  if nargin<1, 
    exp='006_oriZ_dif_110515_DnaNsignal';
    %exp='CM_DnaN-Dif-Gamma-ve-Position1_Series1';
    %exp='CM_HFR1_JK' ; 
    %exp='2a';
  end

  initval=A001_Images_Set_Experiment(user,exp);

  %load the databases--------------------------------------------------
  outname=strcat(initval.basepath,initval.outname{ColourIdx}); %processed inputs
  outname_usr=strcat(initval.basepath,initval.outname_usr);%manual inputs
  if actions.loaddatabase
    load(outname,'S');
    load(outname_usr,'M');
  end
  %------------------------------------------------------------------

  Nbac=Processing_Measure_Database(S)

  [~,chan_no]=size(S) 

  %[chan_no,~]=size(S)%%I THINK THIS WAS AN ERROR
  thisisSizeS = size(S)
  divtimes=[];
  reptimes=[];
  for ch=1:chan_no  %for each channel

    chan_no-ch;
    BF=S(ch).channels.kymo_BF;
    BW=MakeBinaryEdgeImageFrom(BF, initval.BW_Threshold);
    BW=BW';

    if actions.plot %nargin<1,
      figure;
      subplot(1,2,1); P_Color(BF',500,1000, 'grey');;
      subplot(1,2,2); P_Color(1-BW,500,1000, 'grey');
      % [~]=ginput(1); 
      %     demoloadpth='D:\jkerssemakers\My Documents\BN_All_ActiveProjects\BN_All_Programming\BN_All_CommonMatlab\JK_Matlab CleanTools';
      %     demoloadnm='\DemoData_Bacteria';
      %     save(strcat(demoloadpth,demoloadnm),'kymoBF_demo');    
    end


    [~,repno]=size(S(ch).channels.ReplicationCluster);
    %clean the structure if needed
    if isfield(S(ch).channels, 'AutoDivision')
      S(ch).channels = rmfield(S(ch).channels, 'AutoDivision');
    end

    for m=1:repno  %1) For every replicationcluster....
      %First, for convenience, store labels and existence of parent and
      %daughter cells
      ThisDivInfo=Add_FamilyMembers(S,M,ch,m);
      %----------------------------------------------------------------------
      %(if found) create 'left'and 'right' forks that start at termination and extend to the new
      %daugther init points. The new division should lie in the region encompassed by this fork and the subsequent region.
      %this division marks the end of this bacterium.   
      ThisDivInfo=Add_Forks(S,M,ch,ThisDivInfo);
      %Next, find Division time in an area aftr termination, bounded by the cluster positions of the daughter cels. 
      %The first detection of an edge count as division time. If nothing is
      %detected for the duration of the daughtercell's life, or no daughter
      %cells are found in the database, the division time is set as the
      %termination time
      ThisDivInfo=Add_division_time(S,M,ch,ThisDivInfo,BW);

      %Finally, add the division results to the database

      S(ch).channels.AutoDivision(m).family=ThisDivInfo.family;  
      S(ch).channels.AutoDivision(m).Forks=ThisDivInfo.Forks;  
      S(ch).channels.AutoDivision(m).divtype=ThisDivInfo.divtype; 
      S(ch).channels.AutoDivision(m).divisiontime=ThisDivInfo.divisiontime;
      S(ch).channels.AutoDivision(m).accepted=1; % just to be sure field exists
      S(ch).channels.ReplicationCluster(m).Forks=ThisDivInfo.Forks;

      %pcolor(1-BW); colormap grey, shading flat
      title(strcat('exp',exp,'channel' ,num2str(ch)));
      %[~]=ginput(1);
    end

    %_______________________________________________________________________
    % Next, with all division times stored, we re-loop through the database 
    %to find the birth times. 
    %This is simply the division time of the parent.
    %     If no data is stored of the parent, it is the initiation time
    [~,repno]=size(S(ch).channels.ReplicationCluster);
    for m=1:repno  %1) For every replicationcluster....
      ThisDiv=S(ch).channels.AutoDivision(m);
      ThisDivNw=Add_birth_time(S,M,ch,ThisDiv);    
      %Finally, add to database
      S(ch).channels.AutoDivision(m).birthtime=ThisDivNw.birthtime;  
      S(ch).channels.AutoDivision(m).birthtype=ThisDivNw.birthtype;
    end
    %------------------------------------------------------------------------

    for m=1:repno  %1) For every replicationcluster....
      ThisDiv=S(ch).channels.AutoDivision(m);
      ThisDiv=Add_division_edges(S,M,ch,ThisDiv,BW);       
      %Finally, add to database
      S(ch).channels.AutoDivision(m).edges=ThisDiv.edges;
    end


    %plot & save section------------------------------------------------
    %--------------------------------------------------------------------
    %Show the cycle time as difference between birth and division
    for m=1:repno  %1) For every replicationcluster.... 
      strt=S(ch).channels.AutoDivision(m).birthtype;
      stp=S(ch).channels.AutoDivision(m).divtype;
      ok=((strcmp(strt, 'OK') & strcmp(stp, 'OK')) &...
      S(ch).channels.AutoDivision(m).edges.edgesok);  %only the good ones!
      %ok=1;
      if ok
        d=S(ch).channels.AutoDivision(m).divisiontime-S(ch).channels.AutoDivision(m).birthtime;
        r=S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(end)-S(ch).channels.ReplicationCluster(m).PosKyTracCom.frames(1);
        divtimes=[divtimes d];
        reptimes=[reptimes r];
      end
    end



    %Finally, save stuff
    if actions.savedata
      save(outname,'S', '-append');
    end
    display('done');
  end
  [good,bad]=Processing_Measure_GoodDatabase(S)


  figure;
  subplot(1,2,1);
  plot(reptimes,'o'); hold on;
  plot(divtimes,'ro'); hold on;
  %legend(['Replication', 'Automated Division Meas']);
  xlabel('cell');
  ylabel('time, frames');
  subplot(1,2,2);
  plot(reptimes,divtimes, 'k*'); hold on;
  plot([0 110], [0 110], 'k-');
  xlabel('replication time, frames');
  ylabel('division time, frames');
  end


