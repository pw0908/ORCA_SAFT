from ORCA_SAFT import GenerateInput, RunORCA, CleanTemp, CleanOE, GenerateJobFile, ReadPES, ReadSPE, GetPES, GetSPE, MergePES
import os, time
from math import exp
from shutil import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import const

#-----------------------------------------------------------------------

# def RunPure(*, smiles, method, path_orca, mem=62, ncpus=8, orcaproc = None, email = None):
#     # generate input file by first initiallizing a GenerateInput class, with the desired SMILES representation molecule
#     # and methods. The method txt file should contain all the necessary settings. 
#     # pVDZ, pVTZ, pVQZ, pV5Z and pV6Z are all available
#     GI = GenerateInput(smiles, method=method)
    
#     GI.InputPure(mode="parallel",suffix="pure",ncpus=ncpus,orcaproc=orcaproc)
    
#     # run ORCA in parallel mode
#     RunORCA(smiles+"_pure"+"_"+method+".inp",mode="parallel",mem=mem,\
#         ncpus=ncpus,path_orca=path_orca,email=email)

def RunDimer(*, smiles, method, path_orca, rval, sys = "ICLChemEng", hour = 0, \
minute = 30, node = 1, mem = None, ncpus = 8, env = "os", \
orcaproc = None, group = None, ori = None, email = None):
    # generate input file by first initiallizing a GenerateInput class, with the desired SMILES representation molecule
    # and methods. The method txt file should contain all the necessary settings. 
    # pVDZ, pVTZ, pVQZ, pV5Z and pV6Z are all available
    GI = GenerateInput(smiles, method=method)
    
    # get input for pure single molecule
    GI.InputDimer(rval,mode="parallel",suffix="dimer",\
        group=group,ori=ori,orcaproc=orcaproc)
    if None not in {group, ori}:
        RunORCA(smiles+"_dimer"+"_"+method+"_"+group+"_"+ori+".inp",mode="parallel",mem=mem,\
        sys=sys,hour=hour,minute=minute,ncpus=ncpus,\
        env=env,path_orca=path_orca,email=email)
    else:    
        RunORCA(smiles+"_dimer"+"_"+method+".inp",mode="parallel",mem=mem,\
        sys=sys,hour=hour,minute=minute,ncpus=ncpus,\
        env=env,path_orca=path_orca,email=email)

def RunPureCP(*, smiles, method, path_orca, rval, sys = "ICLChemEng", hour = 0, \
minute = 30, node = 1, mem = None, ncpus = 8, env = "os", molecule, \
orcaproc = None, group = None, ori = None, email = None):
    # generate input file by first initiallizing a GenerateInput class, with the desired SMILES representation molecule
    # and methods. The method txt file should contain all the necessary settings. 
    GI = GenerateInput(smiles, method=method)
    
    # get input for pure single molecule
    GI.InputPureCP(rval,molecule=molecule,mode="parallel",suffix="pure",\
        ncpus=ncpus,group=group,ori=ori,orcaproc=orcaproc)
    if None not in {group, ori}:
        RunORCA(smiles+"_pure_"+"CP_"+str(molecule)+"_"+method+"_"+group+"_"+ori+".inp",\
        sys=sys,hour=hour,minute=minute,ncpus=ncpus,\
        env=env,path_orca=path_orca,email=email,mode="parallel",\
        mem=mem)
    else:    
        RunORCA(smiles+"_pure_"+"CP_"+str(molecule)+"_"+method+".inp",\
        sys=sys,hour=hour,minute=minute,ncpus=ncpus,\
        env=env,path_orca=path_orca,email=email,mode="parallel",\
        mem=mem)  

#--------------------------------------------------------------------------
# def PESSPE(*, Nr, smiles, method):
#     ReadPES(Nr = Nr, OutputFile = smiles+"_dimer"+"_"+method+".out")
#     ReadSPE(OutputFile = smiles+"_pure"+"_"+method+".out")

#------------------------------------------------------

# argument functions
def factorial(n):
    fact = 1
    for i in range (1,int(n)+1):
        fact = fact * i
    return fact

# JHBV Ar potential Mol. Phys. 107:20, 2181-2188 (2009), note there's an error in C16. Corrected here.
# length in nm and energy in K
def Ar_JHBV(R):
    A =    4.61330146E7
    a1 =  -2.98337630E1
    # a1 =   2.98337630
    a2 =  -9.71208881
    # a2 =   9.71208881E-2
    a_1 =  2.75206827E-2
    # a_1 =  2.75206827E-1
    a_2 = -1.01489050E-2
    # a_2 = -1.01489050
    b =    4.02517211E1
    # b =    4.02517211
    C =   [4.42812017E-1, # C6
           3.26707684E-2, # C8
           2.45656537E-3, # C10
           1.88246247E-4, # C12
           1.47012192E-5, # C14
           1.17006343E-6] # C16
           # 1.70063432E-6 # C16 paper value
    eps  = 143.12
    Reps = 0.3762
    sig  = 0.3357
    R2 = R*R
    Vexp = A * exp(a1*R + a2*R2 + a_1/R + a_2/R2)
    Vpoly = 0
    for iN in range(3,9):
        Vpoly += C[iN-3]/(R2**iN) * \
        (1 - exp(-b*R) * sum([(b*R)**ik/factorial(ik) for ik in range(0,2*iN+1)]))
    V_Ar = Vexp-Vpoly
    return V_Ar

# Waldrop et al. Kr potential, JCP 142, 204307 (2015) 
# Energy in Hartree and distance in a_0
def Kr_Waldrop(R):
    A =  467.771557
    B =  -43.111875
    C = -509.601417
    alpha = 1.566575
    beta  = 4.083794
    C_array = [
         126.790499, # C6
        5268.109217  # C8
    ]
    A_sh = 1296.0
    alpha_sh = 3.067950
    beta_sh  = 0.3240714
    
    R2 = R*R
    if R<3.40150703:
        V_Kr = (A_sh / R) * exp(-alpha_sh*R + beta_sh*R2)
    else:
        Vexp = (A + B*R + C/R) * exp(-alpha*R)
        Vpoly = 0
        for iN in range(3,5):
            Vpoly += C_array[iN-3]/(R2**iN) * \
            (1 - exp(-beta*R) * sum([(beta*R)**ik/factorial(ik) for ik in range(0,2*iN+1)]))
        V_Kr = Vexp-Vpoly
    return V_Kr

def Mie(x,sigma,eps,r):
    V_Mie = eps*(x[0]/(x[0]-x[1]))*(x[0]/x[1])**(x[1]/(x[0]-x[1]))*((sigma/r)**x[0]-(sigma/r)**x[1])
    return V_Mie

#------------------------------------------------------

email = "tz1416@ic.ac.uk"
smiles = "[Ne]"
group = "Td"
GeometryList = ["AB", "AC", "AD", "CA", "CB", "AA", "EF", "EE", "AE", "AF", "FA", "CF"]
r_0 = 2.9
rint = (2.0-0.9)/49
rval = [2.0*r_0 - rint*r_0*i for i in range(50)]
Nr = 50
path_orca = "../../../../orca_4_2_1_linux_x86-64_openmpi314/orca"
method = "DLPNO-CCSDT_pVQZ_1E-3TCutDO"

folder = "./"+smiles+"_"+method
folder = "./PESSPE"
if not os.path.exists(folder):
    os.mkdir(folder)
os.chdir(folder)
 
copy("../"+method+".txt", "./")
copy("../ORCA_SAFT.py", "./")
copy("../Test.py", "./")
copy("../const.py", "./")
# on cx1 clusters
RunPureCP(smiles=smiles, method=method, path_orca=path_orca, rval=rval, \
minute = 0, node = 1, mem = 90, ncpus = 32, env = "abinitioSAFT", sys = "ICLcx1", hour = 24,molecule = 1, \
orcaproc = 4)
RunDimer(smiles=smiles, method=method, path_orca=path_orca, rval=rval, \
minute = 0, node = 1, mem = 90, ncpus = 32, env = "abinitioSAFT", sys = "ICLcx1", hour = 24, \
orcaproc = 2)
    
#---------------------------------------------------------
# os.chdir("./PESSPE")

# MergePES(smiles+"_pure_CP_1_"+method+"_1-48"+"_PES.pes", smiles+"_pure_CP_1_"+method+"_49"+"_PES.pes", smiles+"_pure_CP_1_"+method+"_PES.pes")
# time.sleep(5)
# MergePES(smiles+"_pure_CP_1_"+method+"_PES.pes", smiles+"_pure_CP_1_"+method+"_50"+"_PES.pes")
# time.sleep(3)
# MergePES(smiles+"_dimer_"+method+"_PES.pes", smiles+"_dimer_"+method+"_17-40_"+"PES.pes")
# time.sleep(3)
# MergePES(smiles+"_dimer_"+method+"_PES.pes", smiles+"_dimer_"+method+"_41-46_"+"PES.pes")
# time.sleep(3)
# MergePES(smiles+"_dimer_"+method+"_PES.pes", smiles+"_dimer_"+method+"_47-50_"+"PES.pes")
# exit()

# ReadPES(Nr = Nr, OutputFile = smiles+"_pure_CP_1_"+method+".out")
# ReadPES(Nr = Nr, OutputFile = smiles+"_dimer_"+method+".out")

# r1, PES1 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
# r3, PES3 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
# for iR in range(len(r3)):
#     PES3[iR] -= 2.0* PES1[iR]
#     PES1[iR] *= const.Hartree2K

# plt.rc('font',family='serif',size=16)
# plt.rc('text',usetex=True)
# plt.rc('xtick',labelsize='small')
# plt.rc('ytick',labelsize='small')
# fig = plt.figure(figsize=(7,5.25))
# ax = fig.add_subplot(1,1,1)
# ax.plot(r3,PES1,color='r',ls='solid',label = "Ne single DLPNO-CCSD(T)/QZ",linewidth=2)
# ax.set_xlabel(r'$r / \AA$')
# ax.set_ylabel(r'$\phi(r)$ / K')
# # plt.ylim((-250,200))
# ax.legend(frameon=False)
# plt.gcf().subplots_adjust(left=0.15)
# plt.savefig(smiles+"_"+method+"_single_CP.png")

#---------------------------------------------------------

# CleanOE("[Ar]")
# CleanTemp("[Ar]")
