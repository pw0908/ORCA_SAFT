from ORCA_SAFT import ORCA_SAFT

# Initialise the ORCA_SAFT module by specifying a smiles and method txt.
# The method txt file should contain all the necessary settings. 
# pVDZ, pVTZ, pVQZ, pV5Z and pV6Z are all available
r = ORCA_SAFT("C",method="HFLD_pVDZ.txt",path_orca="orca")

# test new functions

# Obtain the SAFT parameters from QC calculations
# r.SAFT_Params()

# Plot the resulting potentials (must run SAFT_Params prior to this)
# r.Plotting()

# Clean will remove any undesirable files from the repository
# r.Clean()

# r.MolVol()