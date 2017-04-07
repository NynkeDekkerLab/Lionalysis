E %each cell is an experiment (corresponds to folder)
    E(i,1) %fcelli oufti number or something, Datastruct
        .DataStruct foreach cell, foreach colour channel: 
            x
                amplitude intensity
                x positions
                std x position
                y position
                std y position
                gaussian volume / integrated intensity
                Full cell intensity
            SNR INF because oufti mask
            ydatacrpdR1 Bacpics; see e.g. 
                for i=1:3; figure(i); imagesc( E{1,1}.DataStruct(i,1).ydatacrpdR1{1}); title('1 cfp 2 yfp 3 rfp'); end

            pixels nobody knows Mark's invention
            bx  tracks what spots are rejected by the USER
            Ld x in long-axis short axis  (Does not contain errors, but are the average of l(d) over d and d(l) over l. Can be found by looking at function projectToMesh)            
            CellLength Long-axis in pixels
            Lnorm   is long-axis divided by long-axis maximum in experiment
            BacMesh Oufti Mesh
            
