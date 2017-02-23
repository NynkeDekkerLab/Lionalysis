
function BW=MakeBinaryEdgeImageFrom(KBF,tol);
  % This function transforms an image into an edge-detected, binary image;
  % this is useful if peaks with strongly varying intensities exist. JK13

  KBF=Enhance_FilamentousStructures(KBF,tol);

  BW=0*KBF;
  [r,c]=size(KBF);

  points=[];
  for i=1:r
  prf=KBF(i,:);
  mxi=F110_Get1DpeaksFlatBottom(prf,tol); %get peaks per line
  BW(i,mxi)=1;
  if length(mxi)>0
  ptsi=zeros(length(mxi),3);
  ptsi(:,1)=mxi;
  ptsi(:,2)=i;
  ptsi(:,3)=1;
  points=[points; ptsi];
  end
  %BW(i,2)=1;   %just to have points in every 'frame'
  %BW(i,end-2)=1;
  prf=BW(i,:);
  lm=length(mxi); 
  end


  subplot(1,2,2); P_Color(1-BW,500,1000, 'grey');
  %Series of binary operations
  BW=bwmorph(BW,'clean');      
  BW=bwmorph(BW,'dilate',1);
  %BW=bwmorph(BW,'clean');
  %BW=bwmorph(BW,'erode',1);
  BW=bwmorph(BW,'skel', Inf); 

  %     P_Color(1-BW,500,1000, 'grey');
  % [~]=ginput(1);
end
