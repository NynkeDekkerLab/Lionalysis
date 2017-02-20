
  function filim=Enhance_FilamentousStructures(KBF,tol);
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
       fil_initval.GaussLongAxis=7;
       fil_initval.GaussShortAxis=3;
       fil_initval.GaussSearchRadius=8;
       fil_initval.weighpointsbyintensity=1;     
       filim=JKD2_IM_filament_it(points,KBF,fil_initval);
       filim=filim';
       if 0
        subplot(1,2,1); pcolor(KBF); shading flat; colormap hot
        subplot(1,2,2); pcolor(filim); shading flat; colormap hot
        [~]=ginput(1);
       end
  end
