function BW=MakeBinaryEdgeImagefrom(KBF,tol);
% This function transforms an image into an edge-detected, binary image;
% this is useful if peaks with strongly varying intensities exist. JK13

  BW=0*KBF;
  [r,c]=size(KBF);
  for i=1:r
    prf=KBF(i,:);
    mxi=F110_Get1DpeaksFlatBottom(prf,tol); %get peaks per line
    BW(i,mxi)=1; 
    prf=BW(i,:);
    lm=length(mxi); 
  end
  %Series of binary operations
  BW=bwmorph(BW,'clean');      
  BW=bwmorph(BW,'dilate',1); 
  BW=bwmorph(BW,'skel', Inf); 

  P_Color(1-BW,500,1000, 'grey'); 
end


