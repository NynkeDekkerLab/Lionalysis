
%KymoMaker

%This code  allows the user to click a region with a blob and
%tracks a vector from it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%NOTES: 
%'dipimage' should be installed and started ('dipstart') to make this
%code work
%import an 'images' database first! such as 'Images800.mat'
%Jacob Kerssemakers 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Below, you can define your own paths and files-------------------------

function [initval] = A040_Images_KymoMaker(user,exp)

    if nargin<2
        exp='Exp001_DnaN_TUS_dif_01092016_M';
    end
    if nargin<1
        user = 'MarkPC';
    end

    initval=A001_Images_Set_Experiment(user,exp);

    actions.getdatabase=1;          % default=1 

    for I=1:size(initval.viewchan,2) 
        % Select the opening of the growth channel.
        
        actions.reloadclicks=0; 

        if exist( sprintf('%s%s.mat', initval.basepath, initval.outname_usr))
            action.reloadclicks = 1;
        end
        %Get raw data--------------------------
        if actions.getdatabase 
            ImagesWorkspaceName=strcat(initval.basepath,'Images',num2str(initval.maxfile),'_',initval.outname{I},'.mat');
            fprintf('Loading workspace: %s .\n', ImagesWorkspaceName);
            load(ImagesWorkspaceName,'aa', 'ff');
        end

        % get drift from txt file
        driftpath=strcat(initval.basepath,initval.driftfile);
        drift=dlmread(driftpath);
        %--------------------------------------
         
        [~,~,le]=size(aa); 
        initval.maxfile=le;
        close all;

        im1=squeeze(double(dip_array(aa(:,:,1))));
        if ~initval.correctdrift
            drift=0*drift;
        end
        fr_drift=drift(1,:);

        %Manual click channels - option to re-load
        inname=strcat(initval.basepath,initval.outname_usr);  %separate database storing user actions (limits re-clicking)
        if actions.reloadclicks
        %option al second-channel loading 
            if 1
                load(inname,'manypoints'); 
            else %explicit other click-file!
                inname2=strcat(initval.basepath,'DnaN_TUS_dif_UserInputs');
                load(inname2,'manypoints'); 
            end

        else
            manypoints=Processing_PickManyChannels(im1,initval,'Pick the opening of the desired growth channel.');
            if actions.reloadclicks
                save(inname, 'manypoints', '-append');
            else
                save(inname, 'manypoints');
            end
        end


        %% prepare 'presets': corrected start position, reference map or curves etc.

        for i=1:initval.channelno;
            close all 
            endpoints=[[manypoints(i,1) manypoints(i,2)]; [manypoints(i,3) manypoints(i,4)]];
            presets.twopoints=endpoints;

            presets.type='BF';
            presets.adjustxy=1;
            presets.showmap=0;
            presets.storeref=1;
            presets.useref=0;

            [~,presets,~]=Get_Channel_Map(im1, fr_drift, initval,presets);  %process brightfield; use to re-adjust clicked point

            presets.adjustxy=0;
            presets.storeref=1;

            [~,presets,~]=Get_Channel_Map(im1, fr_drift, initval,presets);  % again process brightfield (show&store adjusted result);

            %% now choose kymograph channel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            presets.type='BF';
            % presets.twopoints=Processing_Pick_Channel(im1,initval,'Choose target channel');

            presets.adjustxy=1;
            presets.storeref=0;
            presets.useref=0;
            [~, presets,~]=Get_Channel_Map(im1, fr_drift, initval,presets);  %process brightfield; re-adjust clicked point

            presets.adjustxy=0;
            presets.useref=1;
            [~]=Get_Channel_Map(im1, fr_drift, initval,presets);   % again process brightfield (show adjusted result);;

            %% Finally, build kymographs using  prior pre-sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            figure;
            presets.showmap=0;
            presets.useref=1;
            [kymo_FL,kymo_BF,chanstk_BF,chanstk_FL]= Build_Kymos(aa,ff,drift,initval,presets);

            close all;
            kymo_BF=fliplr(kymo_BF);
            kymo_FL=fliplr(kymo_FL);
            figure; pcolor(kymo_BF); shading flat; colormap hot;
            figure; pcolor(kymo_FL); shading flat; colormap hot;

            %Savefiles----------------------------------------------------------------
            initval.WorkspaceOutName{i}{I}=strcat(exp,'Chan_x',num2str(ceil(presets.twopoints(1,1))),initval.viewchan{I}); %channel data base
            
            MatOutname = strcat(initval.WorkspaceOutName{i}{I},'.mat');
            lbl1=strcat(initval.basepath,MatOutname);  %path+channel database
            save(lbl1, 'endpoints', 'presets' ,'initval', 'kymo_FL','kymo_BF','chanstk_BF','chanstk_FL');

            %Here we first verify whether we need to create the folder
            FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'Kymographs/',initval.viewchan{I}),'dir');
            if FolderExistence == 0
                mkdir(strcat(initval.basepath,initval.FiguresFolder,'Kymographs/',initval.viewchan{I}))
            end
            if ~ exist(strcat(initval.basepath,'/Kymographs/',initval.viewchan{I}),'dir')
                mkdir(strcat(initval.basepath,'/Kymographs/',initval.viewchan{I}));
            end
            lbl3=strcat(initval.basepath,'/Kymographs/',initval.viewchan{I},'/Kymograph_FL',initval.WorkspaceOutName{i}{I},'.tif'); %kymograph
            lbl4=strcat(initval.basepath,'/Kymographs/',initval.viewchan{I},'/Kymograph_BF',initval.WorkspaceOutName{i}{I},'.tif'); %kymograph
            %For writing Kymos to the Figures folder.
            lbl3_Fig=strcat(initval.basepath,initval.FiguresFolder,'Kymographs/',initval.viewchan{I},'/Kymograph_FL',initval.WorkspaceOutName{i}{I},'.tif'); %kymograph
            lbl4_Fig=strcat(initval.basepath,initval.FiguresFolder,'Kymographs/',initval.viewchan{I},'/Kymograph_BF',initval.WorkspaceOutName{i}{I},'.tif'); %kymograph

            kymim1 = uint8(round(kymo_FL/max(kymo_FL(:))*255 - 1));
            imwrite(kymim1,lbl3,'tif');
            kymim1 = uint8(round(kymo_BF/max(kymo_BF(:))*255 - 1));
            imwrite(kymim1,lbl4,'tif');


            imwrite(kymim1,lbl3_Fig,'tif');
            imwrite(kymim1,lbl4_Fig,'tif');
        end 
        fprintf('Kymos are done for %s.\n',initval.viewchan{I}); 
    end
    fprintf('Finished KymoMaker.\n'); 
end