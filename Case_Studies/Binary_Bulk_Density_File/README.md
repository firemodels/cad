This directory contains the binary input file of vegetation bulk density used in the `fds/Verification/WUI/bulk_density_file.fds` verification case.
The bulk density of canopy foliage used in the example is shown in the image below. Note that only the canopy and not the surface fuel is considered in this case. For verification purposes, the total mass of canopy foliage imported to FDS should be 350 kg.

<img src="example_canopy.png" alt="canopy data for verification" width="400"/>

The script `create_bulk_density_file.py` is used to create the binary file `canopy_foliage.bdf`, following the formatting described in the User Guide. It can easily be modified to proess any 3D array of bulk density data into a binary input file for FDS.