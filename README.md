# Interface Hyper Reduction Study

This repository contains the code that was the basis of a prevous study on reducing the interface mesh of a jointed structure. Please cite the relevant paper: 
```
@article{hyper_paper,
	title = {Reduced order modeling for the dynamics of jointed structures through hyper-reduced interface representation},
	volume = {149},
	issn = {0888-3270},
	doi = {10.1016/j.ymssp.2020.107249},
	journal = {Mechanical Systems and Signal Processing},
	author = {Balaji, N.~N. and Dreher, T. and Krack, M. and Brake, M.~R.~W.},
	month = feb,
	year = {2021},
	pages = {107249},
}
```

This code is a fork of that project code for continued use by the Tribomechadynamics Lab. 
Some additional documentation will be provided here on how to utilize the hyper-reduction code. 
Some folders referenced here may not be publicly available. 


# Steps to Generate Reduced Mesh

1. You will first need to generate a set of matrices describing the system without model reduction. A tutorial on how to do that has been published [here](https://nidish96.github.io/Abaqus4Joints/). This description assumes that coordinates are exported for both sides of the interface. The tutorial may be updated to where that is no longer needed. 
2. Clone this repository. 
3. Under 'FULLJOINT_ROM_PREPARE' create a directory named 'MATS' and download the descriptions of the meshes from the folder [here](https://rice.app.box.com/folder/245530630890). These files should just contain a few items, but it appears that some also contain full M and K matrices that will be ignored for the rest of this process. 
4. TODO: I need to check a different folder on box about having the correct inputs for this step.
4. Open the MATLAB script 'FULLJOINT_ROM_PREPARE/prepare_roms_consint.m'
5. Copy your .mat file produced from Abaqus into the folder './MATRIX_EXTRACTION/RUNS/[SUBFOLDER_NAME]' where '[SUBFOLDER_NAME]' is a subfolder that you create. Copy this subfolder name intothe variable list 'SETDIRS' on line 8 of 'FULLJOINT_ROM_PREPARE/prepare_roms_consint.m'
6. Set the variable 'setid' in this script to select the index of SETDIRS corresponding to your folder.
7. Rename your .mat file to be called 'BRB_WOPRES_MAT.mat' in your folder. 
8. First have to prepare your matrices with a different script, then you can run this one...
9. Need a directory called DATS somewhere and ROMS
10. It looks like 'FULLJOINT_ROM_PREPARE/preparemats.m' was also modified last time I ran this, so that may needed before running the other script.
11. You have to set the method and number of elements in prepare_rom_consint.m. 
12. Set filenames similar to before in preparemats.m
13. New commit here on preparemats.m to change to relative coordinates (may not be needed when doing relative coordinates direct from abaqus). 
