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

# JHBV Ne potential Mol. Phys. 106:1, 133-140 (2008).
# length in nm and energy in K
def Ne_JHBV(R):
    A =    4.02915058383E7
    a1 =  -4.28654039586E1
    a2 =  -3.33818674327
    a_1 = -5.34644860719E-2
    a_2 =  5.01774999419E-3
    b =    4.92438731676E1
    C =   [4.40676750157E-2, # C6
           1.64892507701E-3, # C8
           7.90473640524E-5, # C10
           4.85489170103E-6, # C12
           3.82012334054E-7, # C14
           3.85106552963E-8] # C16
    eps  = 42.152521
    Reps = 0.30894556
    sig  = 0.27612487
    R2 = R*R
    Vexp = A * exp(a1*R + a2*R2 + a_1/R + a_2/R2)
    Vpoly = 0
    for iN in range(3,9):
        Vpoly += C[iN-3]/(R2**iN) * \
        (1 - exp(-b*R) * sum([(b*R)**ik/factorial(ik) for ik in range(0,2*iN+1)]))
    V_Ne = Vexp-Vpoly
    return V_Ne

def Mie(x,sigma,eps,r):
    V_Mie = eps*(x[0]/(x[0]-x[1]))*(x[0]/x[1])**(x[1]/(x[0]-x[1]))*((sigma/r)**x[0]-(sigma/r)**x[1])
    return V_Mie

#--------------------------------------------------------------------
# plot configs

lw = 2.2
labelsize = 16
ticksize = 16
legendsize = 14 #16
rticks = [2.5,3.0,3.5,4.0,4.5,5.0,5.5]
rlabels= ['2.5','3.0','3.5','4.0','4.5','5.0','5.5']

os.chdir('./PESSPE')
smiles = '[Ne]'

#--------------------------------------------------------------------
# plot

method = 'DLPNO-CCSDT_pV5Z'
r7, PES7 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
r8, PES8 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r7)):
    PES8[iR] -= 2.0*PES7[iR]
    PES8[iR] *= const.Hartree2K

method = 'CCSDT_pV5Z_pure'
r9, PES9 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
r10, PES10 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r7)):
    PES10[iR] -= 2.0*PES9[iR]
    PES10[iR] *= const.Hartree2K

method = 'HFLD_pV5Z_conservative'
r11, PES11 = GetPES(smiles+"_pure_CP_1_"+method+"_PES.pes")
r12, PES12 = GetPES(smiles+"_dimer_"+method+"_PES.pes")
for iR in range(len(r7)):
    PES12[iR] -= 2.0*PES11[iR]
    PES12[iR] *= const.Hartree2K

r8m = r8 +[2.5]

plt.rc('font',family='serif')
plt.rc('text',usetex=True)
plt.rc('xtick',labelsize=ticksize)
plt.rc('ytick',labelsize=ticksize)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(r8,PES8,color='r',ls='dotted',label = "DLPNO-CCSD(T)/5Z",linewidth=lw)
ax.plot(r10,PES10,color='g',ls='dashed',label = "CCSD(T)/5Z",linewidth=lw)
ax.plot(r12,PES12,color='b',ls='dashdot',label = "HFLD/5Z",linewidth=lw)
ax.plot(r8,[Mie([9.6977,6.0000],2.8019, 29.875,i) for i in r8[:]],color='m',dashes=[2,2,2,2],label = r"Dufal $\textit{et al.}$",linewidth=lw)
ax.plot(r8,[Mie([13.000,6.0000],2.7760, 37.716,i) for i in r8[:]],color='orange',dashes=[3,1,1,1],label = r"Aasen $\textit{et al.}$",linewidth=lw)
ax.plot(r8m,[Ne_JHBV(i/10) for i in r8m],color='k',ls='solid',label = r"JHBV",linewidth=lw)
ax.set_xlabel(r'$r / \mathrm{\AA}$', fontsize=labelsize)
ax.set_ylabel(r'$\phi(r)$ / K', fontsize=labelsize)
plt.xticks(ticks=rticks,labels=rlabels)
plt.ylim(-50,50)
plt.xlim(2.5,5.5)
ax.legend(frameon=False, fontsize=legendsize, loc='upper right')
plt.gcf().subplots_adjust(left=0.15,bottom=0.13)
# plt.showfig()
plt.savefig("[Ne]_methods_5Z.png")
plt.savefig("[Ne]_methods_5Z.pdf")
plt.savefig("[Ne]_methods_5Z.svg")
