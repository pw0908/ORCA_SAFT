from ORCA_SAFT import GetPES
import os
from math import exp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import const

#--------------------------------------------------------------------

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

def Mie(x,sigma,eps,r):
    V_Mie = eps*(x[0]/(x[0]-x[1]))*(x[0]/x[1])**(x[1]/(x[0]-x[1]))*((sigma/r)**x[0]-(sigma/r)**x[1])
    return V_Mie

#--------------------------------------------------------------------

lw = 2.2
labelsize = 16
ticksize = 16
legendsize = 16
rticks = [3.0,4.0,5.0,6.0,7.0]
rlabels= ['3.0','4.0','5.0','6.0','7.0']

os.chdir('./PESSPE')
smiles = '[Ar]'

method = 'DLPNO-CCSDT_pVDZ'
r1, PES1 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
r2, PES2 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r2)):
    PES2[iR] -= 2.0*PES1[iR]
    PES2[iR] *= const.Hartree2K

method = 'DLPNO-CCSDT_pVTZ'
r3, PES3 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
r4, PES4 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r2)):
    PES4[iR] -= 2.0*PES3[iR]
    PES4[iR] *= const.Hartree2K

method = 'DLPNO-CCSDT_pVQZ'
r5, PES5 = GetPES(smiles+"_pure_CP_1_"+method+"_PES_old.pes")
r6, PES6 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r2)):
    PES6[iR] -= 2.0*PES5[iR]
    PES6[iR] *= const.Hartree2K

method = 'DLPNO-CCSDT_pV5Z'
r7, PES7 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
r8, PES8 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r2)):
    PES8[iR] -= 2.0*PES7[iR]
    PES8[iR] *= const.Hartree2K

plt.rc('font',family='serif')
plt.rc('text',usetex=True)
plt.rc('xtick',labelsize=ticksize)
plt.rc('ytick',labelsize=ticksize)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(r2,PES2,color='r',ls='dotted',label = "DZ",linewidth=lw)
ax.plot(r4,PES4,color='g',ls='dashed',label = "TZ",linewidth=lw)
ax.plot(r6,PES6,color='b',ls='dashdot',label = "QZ",linewidth=lw)
ax.plot(r8,PES8,color='m',dashes=[2,2,2,2],label = "5Z",linewidth=lw)
ax.plot(r2[:],[Ar_JHBV(i/10) for i in r2[:]],color='k',ls='solid',label = r"JHBV",linewidth=lw)
ax.set_xlabel(r'$r / \mathrm{\AA}$', fontsize=labelsize)
ax.set_ylabel(r'$\phi(r)$ / K', fontsize=labelsize)
plt.xticks(ticks=rticks,labels=rlabels)
plt.ylim(-150,150)
plt.xlim(3.0,7.0)
ax.legend(frameon=False, fontsize=legendsize)
plt.gcf().subplots_adjust(left=0.15,bottom=0.13)
# plt.showfig()
plt.savefig("[Ar]_DLPNO_NZ.png")
plt.savefig("[Ar]_DLPNO_NZ.pdf")
plt.savefig("[Ar]_DLPNO_NZ.svg")

#--------------------------------------------------------------------

method = 'CCSDT_pV5Z'
r9, PES9 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
r10, PES10 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r2)):
    PES10[iR] -= 2.0*PES9[iR]
    PES10[iR] *= const.Hartree2K

method = 'HFLD_pV5Z_conservative'
r11, PES11 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
r12, PES12 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r2)):
    PES12[iR] -= 2.0*PES11[iR]
    PES12[iR] *= const.Hartree2K

plt.rc('font',family='serif')
plt.rc('text',usetex=True)
plt.rc('xtick',labelsize=ticksize)
plt.rc('ytick',labelsize=ticksize)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(r8,PES8,color='r',ls='dotted',label = "DLPNO-CCSD(T)/5Z",linewidth=lw)
ax.plot(r10,PES10,color='g',ls='dashed',label = "CCSD(T)/5Z",linewidth=lw)
ax.plot(r12,PES12,color='b',ls='dashdot',label = "HFLD/5Z",linewidth=lw)
ax.plot(r8,[Mie([12.085,6.0000],3.4038, 117.84,i) for i in r8[:]],color='m',dashes=[2,2,2,2],label = r"Dufal $\textit{et al.}$",linewidth=lw)
ax.plot(r8[:],[Ar_JHBV(i/10) for i in r8[:]],color='k',ls='solid',label = r"JHBV",linewidth=lw)
ax.set_xlabel(r'$r / \mathrm{\AA}$', fontsize=labelsize)
ax.set_ylabel(r'$\phi(r)$ / K', fontsize=labelsize)
plt.xticks(ticks=rticks,labels=rlabels)
plt.ylim(-150,150)
plt.xlim(3.0,7.0)
ax.legend(frameon=False, fontsize=legendsize)
plt.gcf().subplots_adjust(left=0.15,bottom=0.13)
# plt.showfig()
plt.savefig("[Ar]_methods_5Z.png")
plt.savefig("[Ar]_methods_5Z.pdf")
plt.savefig("[Ar]_methods_5Z.svg")
