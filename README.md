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
3. Copy your .mat file from Abaqus into the folder './MATRIX_EXTRACTION/RUNS/[SUBFOLDER_NAME]' where '[SUBFOLDER_NAME]' is a subfolder that you create. Copy this subfolder name into the variable list 'SETDIRS' on line 8 of 'FULLJOINT_ROM_PREPARE/prepare_roms_consint.m' Rename your mat file to be called 'BRB_WOPRES_MAT.mat' (WOPRES = Without Prestress). You will also need to copy Elements.dat and Nodes.dat from the Abaqus procedure to this same folder. 
4. Set the variable 'setid' in FULLJOINT_ROM_PREPARE/prepare_roms_consint.m to select the index of SETDIRS corresponding to your folder. 
5. Set the ROM type by setting sel_method to one of {'P', 'U', 'PD'}. These different methods are all described in the journal paper. Then select the number of elements that you want in the reduced mesh. The options for each reduction type are listed in a comment. E.g., 232 is valid for U, but not for P or PD. The variable Nels that the for loop is for needs to be set to this value. 
6. Download the mesh details for the model you are using to a new directory named 'MATS' under 'FULLJOINT_ROM_PREPARE'. For members of the Tribomechadynamics lab, meshes can be found [here](https://rice.app.box.com/folder/245530630890). These files should just contain a few items, but it appears that some also contain full M and K matrices that will be ignored for the rest of this process. 
7. Make a folder under FULLJOINT_ROM_PREPARE named 'ROMS'. Your output model will be put here with the name like 'ROM_U_232ELS.mat' (232 elements with uniform reduction). 
8. The relative coordinate transformation currently needs to be manually defined. Open your abaqus output (Modelmats.mtx) and check if the bottom or top surface outputs nodes first. Look at comments in the relative coordinate transformation to determine which Trel to use (Xt = top, Xb=bot). Abaqus coordinates = Xrel * [Xt - Xb; Xb; Xinternal]. Failing to do so may flip the sign of the contact model and cause other unknown issues if not fully corrected everywhere.
8. Run FULLJOINT_ROM_PREPARE/prepare_roms_consint.m


## Code changes

These are some changes that I have made to the code to get things to run correctly in 2024. These changes are committed along with some file name options for the models that were being produced at the time. 

7. Change line 46 to not include a transpose of R. It should be 'R = sparse(R);' You also need to remove the transpose for Fv on the line before. 
8. The loaded Nodes.dat file contains 3 columns with z coordinates in the last column. However, this code expects only two columns, so that was modified.  
9. Symmetry is no longer enforced after the null space reduction because enforcing symmetry there was creating negative eigenvalues in M. There may be a deeper issue here. 


## Notes

1. If you want to preserve a different number of fixed interface modes than the original Abaqus run, you should be able to change Nrest on line 95 (in the HCBREDUCE function call) to be the number of modes that you with to preserve.
1. I think preparemats.m is just used to reduce to relative coordinates if you do not want to reduce the mesh. It appears that the relative coordinate transform is already applied in prepare_roms_consint.m. The procedure from running this script is probably the same, but note that you need to set the flag 'eliminateXb' to True to get rid of the bottom coordinates and put in relative coordinates. 
3. To run the matrices resuling from this process in [TMDSimPy](https://github.com/tmd-lab/tmdsimpy), an additional post processing step to save for python is needed. See [here](https://github.com/tmd-lab/microslip-rough-contact/blob/main/EPMC_SIMS/save_matrices_for_py.m) for a script.


## Draft of Steps to do next.


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
