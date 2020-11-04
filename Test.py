from ORCA_SAFT import GenerateInput, RunORCA, CleanTemp, CleanOE, GenerateJobFile, ReadPES, ReadSPE, GetPES, GetSPE, MergePES
import os, time
from math import exp
from shutil import copy
import matplotlib
matplotlib.use('Agg')
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

def RunDimer(*, smiles, method, path_orca, rval, mem=None, ncpus=8, orcaproc = None, group = None, ori = None, email = None):
    # generate input file by first initiallizing a GenerateInput class, with the desired SMILES representation molecule
    # and methods. The method txt file should contain all the necessary settings. 
    # pVDZ, pVTZ, pVQZ, pV5Z and pV6Z are all available
    GI = GenerateInput(smiles, method=method)
    
    # get input for pure single molecule
    GI.InputDimer(rval,mode="parallel",suffix="dimer",ncpus=ncpus,\
        group=group,ori=ori,orcaproc=orcaproc)
    
    RunORCA(smiles+"_dimer"+"_"+method+".inp",mode="parallel",mem=mem,\
        ncpus=ncpus,path_orca=path_orca,email=email)

def RunPureCP(*, smiles, method, path_orca, rval, molecule, mem=None, ncpus=8, orcaproc = None, group = None, ori = None, email = None):
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
smiles = "[Ar]"
r_0 = 3.5
rint = (2.0-0.9)/49
rval = [2.0*r_0 - rint*r_0*i for i in range(50)]
Nr = 50
path_orca = "../../../../orca_4_2_1_linux_x86-64_openmpi314/orca"
method = "CCSDT_pVQZ_F12"
# method = "HFLD_pV5Z_conservative"
# rval = 0.87*r_0
# rval = rval[48:49]
folder = "./"+smiles+"_"+method

# folder += "_SPE2"
# rval = rval[49:50]

# if not os.path.exists(folder):
#     os.mkdir(folder)
# os.chdir(folder)

# copy("../"+method+".txt", "./")
# copy("../ORCA_SAFT.py", "./")
# copy("../Test.py", "./")
# copy("../const.py", "./")

# RunPureCP(smiles=smiles,method=method,path_orca=path_orca,rval=rval,molecule=1,mem=10,ncpus=8,email=email)
# time.sleep(3)
# RunDimer (smiles=smiles,method=method,path_orca=path_orca,rval=rval,           mem=10,ncpus=8,email=email)

# PESSPE(Nr = Nr, smiles = smiles, method = method)
# Nr = 1
# os.chdir("./PESSPE")
# ReadPES(Nr = Nr, OutputFile = smiles+"_dimer_"+method+".out")
# ReadPES(Nr = Nr, OutputFile = smiles+"_pure_CP_1_"+method+".out")


#---------------------------------------------------------
os.chdir("./PESSPE")

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

r1, PES1 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
r2, PES2 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r2)):
    PES2[iR] -= 2*PES1[iR]
    PES2[iR] *= const.Hartree2K

# r3, PES3 = GetPES(smiles+"_pure_CP_1_"+"DLPNO-CCSDT_pV5Z"+"_PES.pes")
# r4, PES4 = GetPES(smiles+"_dimer_"+"DLPNO-CCSDT_pV5Z"+"_PES.pes")
# for iR in range(len(r2)):
#     PES4[iR] -= 2*PES3[iR]
#     PES4[iR] *= const.Hartree2K

# rlist = []
# PESlist = []
# PESlabel = []

# r3, PES3 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
# SPE1 = GetSPE(smiles+"_pure_"+method+"_SPE.spe")

# PES3 = [(i - 2*SPE1) * const.Hartree2K for i in PES3]
# PES4 = [(i - 2*SPE4) * c.Hartree2K for i in PES4]
# PES5 = [(i - 2*SPE5) * c.Hartree2K for i in PES5]
# PES6 = [(i - 2*SPE6) * c.Hartree2K for i in PES6]
# PES7 = [(i - 2*SPE1) * c.Hartree2K for i in PES7]
# print(PES7)
plt.rc('font',family='serif',size=16)
plt.rc('text',usetex=True)
plt.rc('xtick',labelsize='small')
plt.rc('ytick',labelsize='small')
fig = plt.figure(figsize=(7,5.25))
ax = fig.add_subplot(1,1,1)
ax.plot(r2,PES2,color='r',ls='solid',label = "DLPNO-CCSD(T)/5Z with CP",linewidth=2)
# ax.plot(r4,PES4,color='g',ls='solid',label = "DLPNO-CCSD(T)",linewidth=2)
# ax.plot(rval[:],[Mie([12.085,6.0000],3.4038, 117.84,i) for i in rval[:]],color='g',ls='solid',label = "Dufal et al.",linewidth=2)
# ax.plot(r4,PES4,color='g',ls='solid',label = "DLPNO-CCSD(T)/5Z",linewidth=2)
# ax.plot(r5,PES5,color='b',ls='solid',label = "HFLD/6Z",linewidth=2)
# ax.plot(r6,PES6,color='m',ls='solid',label = "CCSD(T)-F12/TZ",linewidth=2)
# ax.plot(r7,PES7,color='orange',ls='solid',label = "CCSD(T)/5Z",linewidth=2)
# ax.plot(rval[:],[Ar_JHBV(i/10) for i in rval[:]],color='k',ls='--',label = "Jager et al.",linewidth=2)
ax.set_xlabel(r'$r / \AA$')
# # ax.set_ylabel(r'$\phi(r)$ / kcal/mol')
ax.set_ylabel(r'$\phi(r)$ / K')
# ax.legend(frameon=False)
plt.gcf().subplots_adjust(left=0.15)
plt.savefig(smiles+"_"+method+".png")

#---------------------------------------------------------

# CleanOE("[Ar]")
# CleanTemp("[Ar]")
