from ORCA_SAFT import GenerateInput, RunORCA, CleanTemp, CleanOE, GenerateJobFile, ReadPES, ReadSPE, GetPES, GetSPE
import os
from math import exp
from shutil import copy
import matplotlib.pyplot as plt
import const

#-----------------------------------------------------------------------

def RunPure(*, smiles, method, path_orca, mem=62, ncpus=8, orcaproc = None, email = None):
    # generate input file by first initiallizing a GenerateInput class, with the desired SMILES representation molecule
    # and methods. The method txt file should contain all the necessary settings. 
    # pVDZ, pVTZ, pVQZ, pV5Z and pV6Z are all available
    GI = GenerateInput(smiles, method=method)
    
    GI.InputPure(mode="parallel",suffix="pure",ncpus=ncpus,orcaproc=orcaproc)
    
    # run ORCA in parallel mode
    RunORCA(smiles+"_pure"+"_"+method+".inp",mode="parallel",mem=mem,\
        ncpus=ncpus,path_orca=path_orca,email=email)

def RunDimer(*, smiles, method, path_orca, rval, mem=62, ncpus=8, orcaproc = None, group = None, ori = None, email = None):
    # generate input file by first initiallizing a GenerateInput class, with the desired SMILES representation molecule
    # and methods. The method txt file should contain all the necessary settings. 
    # pVDZ, pVTZ, pVQZ, pV5Z and pV6Z are all available
    GI = GenerateInput(smiles, method=method)
    
    # get input for pure single molecule
    GI.InputDimer(rval,mode="parallel",suffix="dimer",ncpus=ncpus,\
        group=group,ori=ori,orcaproc=orcaproc)
    
    RunORCA(smiles+"_dimer"+"_"+method+".inp",mode="parallel",mem=mem,\
        ncpus=ncpus,path_orca=path_orca,email=email)

def RunPureCP(*, smiles, method, path_orca, rval, molecule, mem=62, ncpus=8, orcaproc = None, group = None, ori = None, email = None):
    # generate input file by first initiallizing a GenerateInput class, with the desired SMILES representation molecule
    # and methods. The method txt file should contain all the necessary settings. 
    # pVDZ, pVTZ, pVQZ, pV5Z and pV6Z are all available
    GI = GenerateInput(smiles, method=method)
    
    # get input for pure single molecule
    GI.InputPureCP(rval,molecule=molecule,mode="parallel",suffix="pure",\
        ncpus=ncpus,group=group,ori=ori,orcaproc=orcaproc)
    
    RunORCA(smiles+"_pure_"+"CP_"+str(molecule)+"_"+method+".inp",\
        mode="parallel",mem=mem,ncpus=ncpus,path_orca=path_orca,\
        email=email)     

#--------------------------------------------------------------------------
def PESSPE(*, Nr, smiles, method):
    ReadPES(Nr = Nr, OutputFile = smiles+"_dimer"+"_"+method+".out")
    ReadSPE(OutputFile = smiles+"_pure"+"_"+method+".out")

#------------------------------------------------------

email = "tz1416@ic.ac.uk"
smiles = "[Kr]"
r_0 = 3.5
rint = (2.0-0.9)/49
rval = [2.0*r_0 - rint*r_0*i for i in range(50)]
Nr = 50
path_orca = "../../../../orca_4_2_1_linux_x86-64_openmpi314/orca"
method = "DLPNO-CCSDT_pVQZ"
# method = "HFLD_pV5Z_conservative"
# rval = 0.87*r_0
# rval = rval[46:]
folder = "./"+smiles+"_"+method

if not os.path.exists(folder):
    os.mkdir(folder)
os.chdir(folder)

copy("../"+method+".txt", "./")
copy("../ORCA_SAFT.py", "./")
copy("../Test.py", "./")
copy("../const.py", "./")

RunPureCP(smiles=smiles,method=method,path_orca=path_orca,rval=rval,molecule=1,mem=10,ncpus=8,email=email)

# RunPure(smiles = smiles, method = method, path_orca = path_orca, rval = rval, mem=62, ncpus=1, node=1, hour = 0, minute = 30)
# RunDimer(smiles = smiles, method = method, path_orca = path_orca, rval = rval, mem=62, ncpus=32, node=1, hour = 24, minute = 0, orcaproc = 3)
# PESSPE(Nr = Nr, smiles = smiles, method = method)
# Nr = 12
# ReadPES(Nr = Nr, OutputFile = smiles+"_dimer"+"_"+method+"_1-12"+".out")
# ReadSPE(OutputFile = smiles+"_pure"+"_"+method+".out")

#---------------------------------------------------------

# r1, PES1 = GetPES(smiles+"_dimer_"+method+"_1-12"+"_PES.pes")
# r2, PES2 = GetPES(smiles+"_dimer_"+method+"_13-15"+"_PES.pes")
# r3, PES3 = GetPES(smiles+"_dimer_"+method+"_16"+"_PES.pes")
# r4, PES4 = GetPES(smiles+"_dimer_"+method+"_17-40"+"_PES.pes")
# r5, PES5 = GetPES(smiles+"_dimer_"+method+"_41-46"+"_PES.pes")
# r6, PES6 = GetPES(smiles+"_dimer_"+method+"_47-50"+"_PES.pes")

# PES7 = PES1+PES2+PES3+PES4+PES5+PES6
# r7 = r1+r2+r3+r4+r5+r6

# SPE1 = GetSPE(smiles+"_pure_"+method+"_SPE.spe")

# PES1 = [(i - 2*SPE1) * c.Hartree2K for i in PES1]
# PES4 = [(i - 2*SPE4) * c.Hartree2K for i in PES4]
# PES5 = [(i - 2*SPE5) * c.Hartree2K for i in PES5]
# PES6 = [(i - 2*SPE6) * c.Hartree2K for i in PES6]
# PES7 = [(i - 2*SPE1) * c.Hartree2K for i in PES7]
# print(PES7)
# plt.rc('font',family='serif',size=16)
# plt.rc('text',usetex=True)
# plt.rc('xtick',labelsize='small')
# plt.rc('ytick',labelsize='small')
# fig = plt.figure(figsize=(7,5.25))
# ax = fig.add_subplot(1,1,1)
# ax.plot(r1,PES1,color='r',ls='solid',label = "DLPNO-CCSD(T)/QZ",linewidth=2)
# ax.plot(r4,PES4,color='g',ls='solid',label = "DLPNO-CCSD(T)/5Z",linewidth=2)
# ax.plot(r5,PES5,color='b',ls='solid',label = "HFLD/6Z",linewidth=2)
# ax.plot(r6,PES6,color='m',ls='solid',label = "CCSD(T)-F12/TZ",linewidth=2)
# ax.plot(r7,PES7,color='orange',ls='solid',label = "CCSD(T)/5Z",linewidth=2)
# ax.plot(r1,[Ar_JHBV(i/10) for i in r1],color='k',ls='--',label = "CCSD(T)/6Z(JHBV fit)",linewidth=2)
# ax.set_xlabel(r'$r / \AA$')
# # ax.set_ylabel(r'$\phi(r)$ / kcal/mol')
# ax.set_ylabel(r'$\phi(r)$ / K')
# ax.legend(frameon=False)
# # plt.tight_layout()
# plt.gcf().subplots_adjust(left=0.15)
# plt.savefig(smiles+"_"+"CCSD(T)_5Z_Kr"+".png")

#---------------------------------------------------------

# CleanOE("[Ar]")
# CleanTemp("[Ar]")
