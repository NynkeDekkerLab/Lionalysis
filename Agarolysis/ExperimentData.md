E %each cell is an experiment (corresponds to folder)
    E(i,1) %fcelli oufti number or something, Datastruct
        .DataStruct foreach cell, foreach colour channel: 
            x best spot
                amplitude intensity
                x positions
                std x position
                y position
                std y position
                gaussian volume / integrated intensity
                Full cell intensity
            bx all spots
            ld all spots in cell coordinates
            SNR INF because oufti mask
            ydatacrpdR1 Bacpics; see e.g. 
                for i=1:3; figure(i); imagesc( E{1,1}.DataStruct(i,1).ydatacrpdR1{1}); title('1 cfp 2 yfp 3 rfp'); end

            pixels nobody knows Mark's invention 
            CellLength Long-axis in pixels
            Lnorm   is long-axis divided by long-axis maximum in experiment
            BacMesh Oufti Mesh
            
