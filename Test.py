from ORCA_SAFT import GenerateInput, RunORCA, CleanTemp, CleanOE, GeneratePBS, ReadPES

# smiles = "[Ar]"
# r_0 = 3.5
# rint = (2-0.95)/49
# rval = [2*r_0 - rint*r_0*i for i in range(50)]
Nr = 50

# # generate input file by first initiallizing a GenerateInput class, with the desired SMILES representation molecule
# # and methods. The method txt file should contain all the necessary settings. 
# # pVDZ, pVTZ, pVQZ, pV5Z and pV6Z are all available
# GI = GenerateInput(smiles, method="HFLD_pVDZ.txt")

# # get input for pure single molecule
# GI.InputDimer(rval, mode = "parallel", suffix = "dimer", ncpus = 8, group = None, ori = None)

# # run ORCA in parallel mode
# RunORCA(smiles+"_dimer.inp", mode = "parallel", hour = 0, minute = 30, node = 1, mem = 16, ncpus = 8,
# env = "abinitioSAFT", path_orca = "../../../orca_4_2_1_linux_x86-64_openmpi314/orca")

# ReadPES(Nr = Nr, OutputFile = "[Ar]_dimer.out", EnergyType = "MDCI energy")
CleanOE("[Ar]")
CleanTemp("[Ar]")
