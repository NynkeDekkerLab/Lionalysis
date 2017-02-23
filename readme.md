NynkeDekkerLab code for fluorescent Analysis. These are related to our in-vivo fluorescence experiments using agarose pads and microfluidic devices. 

The original contributors to this software framework are:
-  Dr. Jacob Kerssemakers
-  Roy de Leeuw, M. Sc.
-  Mark (S.) Hu, B. Sc.

The code is currently maintained by J. S. de Boer, M. Sc. of the Nynke Dekker Lab.


Required:
- DIPImage and added to path ~DipInstallationDir\dipstart.m
- Run dipstart
- UTrack and added to path ~UTrackInstallationDir\software\ and ~UTrackInstallationDir\software\mex
- NanSuite and add to path ~NanSuiteInstallationDir
- A clone of https://github.com/bmelinden/vbSPT.git added to path

============================
- KymoCode: A module for analysing fluorescence images of bacteria that grow linearly, constraint in a straight growth channel. It consists of several parts that operate in sequence to produce a series of image stacks of a single cell during a period of cell growth in between two divisions, the so called Backpics. It first does image preparation, mainly rolling background subtraction. It hten stores the images in matlab format for easier loading. Next, it calculated drift by correlating a user selected region. Next, it makes the Kymographs. Finally, the user identifies replication events in the Kymographs and Bacpics are created.

- LionFit: BacPics are self-constained single cell images, for which the fluorescent spot analysis software can be seperate. LionFit finds the intensity, position and the uncertainty in the position of the fluorescent spot within the BacPics, as well as the illumination intensity of the entire cell. Uses http://www.nature.com/nmeth/journal/v7/n5/abs/nmeth.1449.html

- Agarolysis: Similarly to the KymoCode, the Agrolysis analyses fluorescence images of bacteria by producing BacPics. However, the experiments are performed on agarose pads instead of in a contrained linear channel. agarolysis also uses LionFit to analyze several flureoscence properties of the bacterial cells. Hwoever, in Agarolysis some post-processing on the LionFit results is needed, as LionFit returns the Cartesian but not the cellular coordinates of hte Fluorescence spots. After setting the configuration by a neat user UI, it calls Oufti with standardised settings (i.e. set in experiment config) to detect the cell outlines in the phase contrast image. It then uses a rolling background subtraction on the images. Next, it makes pictures around the Oufti-detected PC outlines. It then uses LionFit to detect spots in each Bacpics, in principle via a global treshold scheme. Then, again using Oufti, it projects the LionFit coordinates onto the Oufti meshes. Finally, each Bacpic is shown and the user can validate the final results.  It then outputs the found results to the matlab workspace and in Results.mat.

- TigerCreate: Selects the diffusion constant in an user interface created and embedded into BlurLba. In this interface the user can select up to three different diffusion constants, indicate the percentage of spots to behave according to these diffusion constant s and how the distribution changes over the course of the simulation. 

- LionDiffusion: Analyses diffusion properties of fluorescent spots using single-particle tracking using homebrew software, utrack and vbSPTgui.

- vbSPTgui: Uses utrack's trajectories to analyse diffusion properties using a hidden Markov model. Requires each experiment to be in a directory named with the experiment number.

- LionSMC is effectively an administrative tool that converts UTrack results in a histogram that gives information egarding the intensity of a spot, so that you can assess the number of spots involved in a particular result. To use it, first use movieSelectorGui and U-track to get results, then analyse them using LionStat and LionSMC.


=============================
Hints
- Don't go to negative peak tresholds
- If the code crashes, try removing the Edited Images and Backpic folders. They somehow are able to interfere?