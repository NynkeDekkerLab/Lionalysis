 function RoicotrascaKymo(initval,Imgname,Beampath, id)
    %% Presets

    OSslash = initval.OSslash;
    jarpath = initval.Supportpath;
    Imgspath = strcat(initval.basepath,initval.FLpath);
    Rval = initval.flresize;
    Tval = initval.fltrans;
    Raw = initval.rawfolder;

    Outpath = strcat(Imgspath,Imgname,OSslash);

    if ~exist(Outpath,'dir')
        mkdir(Outpath);
    end

    % Add ImageJ and MIJI to path
    javaaddpath(strcat(jarpath,'mij.jar'))
    javaaddpath(strcat(jarpath,'ij.jar'))

    %% Rolling ball
    Rballedname = strcat('Rballed_',Imgname,'.tif');
    % copyfile(strcat(Imgspath,Imgname),strcat(Outpath,Rballedname));
    % disp('Copying file...') 

    % Start the MIJi.m interface between matlab and imageJ (opens a matlab-owned imageJ instance)
    fprintf('Opening image (%d)...\n', id);
    MIJ.start
    % Open imageJ, with base directory set to $open
    MIJ.run('Image Sequence...',strcat('open=[',Imgspath,Imgname,OSslash,Raw,']'))
    % Open a Tiff in imageJ, path is $open/$path
    MIJ.run('Tiff...',strcat('path=[',Outpath,Imgname,'.tif]'))

    fprintf('Performing rolling ball background subtraction for image %d.\n', id);
    % Run the subtract background module, with $rolling='10 stack' indicating a ten point rolling ball background subtraction.
    MIJ.run('Subtract Background...', 'rolling=10 stack')
    % Save the result
    MIJ.run('Tiff...',strcat('path=[',Outpath,Rballedname,']'))
    % Literally what it says
    MIJ.closeAllWindows
    MIJ.exit

    fprintf('Rolling ball finished for image %d\n', id);

    %%

    imgpath = strcat(Outpath,Rballedname);
    imginfo = imfinfo(imgpath);
    num_images = numel(imginfo);
    RITpath = strcat(Outpath,'RIT_',Imgname,'.tif');

    Beamimg = im2double(imread(Beampath));
    maxVAlue= max(max(Beamimg));
    NewBeamImg = Beamimg./maxVAlue;

    RYval = floor(Rval*imginfo(1).Width);
    RXval = floor(Rval*imginfo(1).Height);

    fprintf('Performing illumination correction and transformation for image %d .\n', id);
    for k = 1:num_images;
        fprintf('\t Image %d out of %d. \n', k, num_images);
        img = im2double(imread(imgpath,k,'Info',imginfo));
        
        % Illumunation Correction
        RIimg = double(imdivide(img,NewBeamImg));
        
        %% Transformation
        
        RITimg = RIimg;
        if ~isequal(Rval,1)
            RITimg = imresize(RITimg,[RXval,RYval],'bilinear');
        end
        
        if ~isequal(Tval,[0,0])
            RITimg = imtranslate(RITimg,Tval,'FillValues',0);
        end  
        
        % Write to file
        imwrite(uint16(65535*RITimg),RITpath,'WriteMode','append','Compression','none');
    end

    fprintf('Roicotrasca done for image %d\n', id)
end

