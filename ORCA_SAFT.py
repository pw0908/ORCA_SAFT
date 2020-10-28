import subprocess,os,glob
import numpy as np
from math import exp, sqrt, acos, cos, sin
from scipy.interpolate import interp1d
from scipy.optimize import fsolve,minimize
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
from openbabel import pybel
from copy import deepcopy

# Generate input files for ORCA calculations
# This input file generator is expected to generate input files
# for different purposes, e.g. single molecule, dimer with different 
# orientations and else 
class GenerateInput(object):
 def __init__(self,smiles,method="HFLD_pVDZ",path_multiwfn="multiwfn"):
    
    # Create a Molecule object in pybel
    mol = pybel.readstring("smi",smiles)

    # Set it up in 3D coordinates
    mol.make3D()

    self.smiles = smiles
    self.method = method
    self.path_multiwfn = path_multiwfn
    self.AtomInfo = []
    self.AtomInfo.append(deepcopy(AtomInfo(mol)))
    self.AtomInfo.append(deepcopy(AtomInfo(mol)))
 
 # generate calculation for single molecule
 def InputPure(self, mode = "single", suffix = "pure", ncpus = 8, orcaproc = None):
    if orcaproc is None:
        nprocs = ncpus
    else:
        nprocs = orcaproc
    AI = self.AtomInfo
    PureFileName = self.smiles+"_"+suffix+"_"+self.method+".inp"
    if os.path.exists(PureFileName):
        os.remove(PureFileName)
    file = open(PureFileName,"w")
    with open(self.method+".txt") as f:
        lines = f.read().split("\n")
        for i in range(len(lines)):
            file.write(lines[i]+"\n")
    
    # processor number
    if mode == "parallel":
        file.write("%pal nprocs "+str(int(nprocs))+" end\n")

    file.write("\n")
    
    # coordinates from AtomInfo
    file.write("* xyz 0 1\n")
    for iAtom in range(AI[0].NAtoms):
        file.write("  "+AI[0].type[iAtom]+"(1)\t"+str(AI[0].coord[iAtom][0])+"\t"+str(AI[0].coord[iAtom][1])+"\t"+str(AI[0].coord[iAtom][2])+"\n")
    file.write("*\n\n\n")    
    file.close()

 def InputDimer(self, rval, mode = "single", suffix = "dimer", ncpus = 8, group = None, ori = None, orcaproc = None):
    if orcaproc is None:
        nprocs = ncpus
    else:
        nprocs = orcaproc
    AI = self.AtomInfo
    # so far considered only Td geometry 
    if group == "Td":
        AI[0].TdOrient(ori[0])
        AI[1].TdOrient(ori[1])
    # possibly write a function in future to systematically check the two-body orientations in a 
    # systematic way.
    else:
        pass
    DimerFileName = self.smiles+"_"+suffix+"_"+self.method
    if None not in {group, ori}:
        DimerFileName += ("_"+group+ori)
    DimerFileName += ".inp"
    if os.path.exists(DimerFileName):
        os.remove(DimerFileName)
    file = open(DimerFileName,"w")
    # preset methods
    with open(self.method+".txt") as f:
        lines = f.read().split("\n")
        for i in range(len(lines)):
            file.write(lines[i]+"\n")

    if mode == "parallel":
        file.write("%pal nprocs "+str(int(nprocs))+" end\n")
    
    # coords block
    file.write("\n")
    file.write("%coords\n")
    file.write("CTyp\txyz\n")
    file.write("Charge\t0\n")
    file.write("Mult\t1\n")
    # define parameters for varied r
    file.write("pardef\n")
    if isinstance(rval, float):
        file.write("\tr = "+str(rval)+";\n")
    elif isinstance(rval, list):
        Crval = "["
        for irval in rval:
            Crval += str(irval) + " " 
        Crval += "]" 
        file.write("\tr "+Crval+";\n")
    file.write("end\n")
    file.write("coords\n")
    
    # dimer position
    # Molecule 1
    for iAtom in range(AI[0].NAtoms):
        file.write("  "+AI[0].type[iAtom]+"(1)\t"\
        +str(AI[0].coord[iAtom][0])+"\t"\
        +str(AI[0].coord[iAtom][1])+"\t"\
        +str(AI[0].coord[iAtom][2])+"\n")
        
    # Molecule 2
    for iAtom in range(AI[1].NAtoms):
        file.write("  "+AI[1].type[iAtom]+"(2)\t"\
        +str(AI[1].coord[iAtom][0])+"\t"\
        +"{"+ str(AI[1].coord[iAtom][1]) +"+r}\t"\
        +str(AI[1].coord[iAtom][2])+"\n")

    file.write("end\n")
    file.write("end\n")
    # use the MOs generated in the last iteration
    if isinstance(rval, list):
        file.write("%method\n")
        file.write("ScanGuess MORead\n")
        file.write("end\n")
    file.write("\n\n\n")    
    file.close()

# collection of functions to get SAFT parameters out, and 
# sigma, epsilon -> Vol and m -> lambdas 
# write this as a class? or a separate python file?
class GetSAFTParameters(object):
 def __init__(self):
    pass

# function for running ORCA, default is single mode
# for the meaning of the kwargs, see GeneratePBS() function
def RunORCA(InputFile, mode = "single", \
 hour = 0, minute = 30, node = 1, mem = 16, ncpus = 8, env = "abinitioSAFT", path_orca = "orca"):
    if mode == "single":
        # just run orca
        print(path_orca+" "+InputFile+">"+InputFile.split(".")[0]+".out")
        subprocess.run(path_orca+" "+InputFile+">"+InputFile.split(".")[0]+".out", shell=True)
    if mode == "parallel":
        # check if the path typed in is the full orca path
        if path_orca == "orca":
            raise Exception("Exact orca path need to be specified")
        # generate PBS file, and get PBS file name
        PBSName = GeneratePBS(InputFile=InputFile,hour=hour,minute=minute,node=node,mem=mem,ncpus=ncpus,env=env,path_orca=path_orca)
        print("qsub "+'"'+PBSName+'"')
        subprocess.run("qsub "+'"'+PBSName+'"', shell=True)

# function for cleaning unnecessary temporary ORCA files
def CleanTemp(smiles):
    matches = [".pbs", ".o", ".e", ".out", ".inp", ".pes", ".spe"]
    for filename in glob.glob("./*"):
        # only clean files with the specified SMILES, prevent the case that specifies "C", but clean files such as "CC"
        if smiles in filename:
            if smiles in {filename.split("_")[0][2:], filename.split("_")[1]}:
        	    # Don't remove pbs, potential energy surface, single point energy, HPC output, error, ORCA input and output files
                if not any(CFormat in filename for CFormat in matches):
                    os.remove(filename)

# function for cleaning unnecessary output and error message from HPC
def CleanOE(smiles):
    # Remove files that aren't necessary
    for filename in glob.glob("./*"):
        # only clean files with the specified SMILES, prevent the case that specifies "C", but clean files such as "CC"
        if smiles in filename:
            if smiles in {filename.split("_")[0][2:], filename.split("_")[1]}:
                matches = [".o", ".e"]
        	    # remove HPC output and error files
                if any(CFormat in filename for CFormat in matches):
                    os.remove(filename)

# function for generating PBS file for submitting HPC job
# hour, minute: HPC machine time
# node, mem, ncpus: # of node, memory size (in GB) and # of cpus (also equal to processors)
# env: anaconda3 environment for running ORCA
# path_orca: orca path, need to be the full path, but relative path (../ and ./) is acceptable
def GeneratePBS(*, InputFile, hour, minute, node, mem, ncpus, env, path_orca):
    # extract PBS name from the input ORCA file
    PBSName = "ORCA_"+InputFile.split(".")[0]+".pbs"
    if os.path.exists(PBSName):
        os.remove(PBSName)
    file = open(PBSName,"w")
    # wall time need to be in HH:MM:SS
    if hour < 10:
        Chour = "0"+str(int(hour))
    elif hour > 10:
        Chour = str(int(hour))
    if minute < 10:
        Cminute = "0"+str(int(minute))
    elif minute > 10:
        Cminute = str(int(minute))
    file.write("#PBS -l walltime="+Chour+":"+Cminute+":00\n")
    # resource request
    file.write("#PBS -l select="+str(int(node))+":ncpus="+str(int(ncpus))+":mem="+str(int(mem))+"gb:mpiprocs="+str(int(ncpus))+"\n")  
    file.write("\n")       
    # activate anaconda3 and the environment
    file.write("module load anaconda3/personal\n")
    file.write("module load orca/4.2.1\n")
    file.write("source activate "+env+"\n")
    # need to set an environment variable of OpenMPI
    file.write("export OMPI_MCA_btl=self,vader,tcp\n")
    file.write("\n")
    # the job submitted will be sent to another location to be run, so if we want to access to
    # some locations in home directory, cd to $PBS_O_WORKDIR is needed
    file.write('cd "$PBS_O_WORKDIR" \n')
    # call ORCA to run the input script
    file.write('"'+path_orca+'"'+" "+'"'+InputFile+'"'+">"+'"'+InputFile.split(".")[0]+".out"+'"'+"\n")
    file.write("\n")
    # deactivate the environment in the end
    file.write("conda deactivate\n")
    file.close()
    return PBSName

# read the potential energy surface (PES) data from the output file as specified in OutputFile
# Nr: the # of r calculated
# EnergyType: in some calculations, there are different types of energies, for HFLD, there are
# "'Actual Energy'", "SCF Energy" and "MDCI Energy". The "'Actual Energy'" is the same as 
# "SCF Energy", but they are not the energy we are looking for. We need "MDCI Energy" because
# the CI part between the two fragments are added to the SCF energy. 
def ReadPES(*, Nr, OutputFile):
    # generate PES file name, in .pes
    PESFile = OutputFile.split(".out")[0]+"_PES.pes"
    # a list of PES values
    rList = []
    PESList = []
    iPES = 0
    ir = 0
    with open(OutputFile,"r") as f:
        lines = f.read().split("\n")
        for iLine,line in enumerate(lines):
            # find the correct line for PES data
            if "FINAL SINGLE POINT ENERGY" in line:
                words = line.split()
                if len(words)>6:
                    PESList.append(float(words[6]))
                else:
                    PESList.append(float(words[4]))
                iPES += 1
            if "R  :" in line:
                words = line.split("R  :")
                rList.append(float(words[1]))
                ir += 1
            # else:
            #     if "r =" in line:
            #         words = line.split("r =")
            #         rList.append(float(words[1].split(";")[0]))
            #         ir += 1
    if (iPES != Nr):
        raise Exception("Either the # of r or the # of energy does not match the desired #, please check outputs")

    if os.path.exists(PESFile):
        os.remove(PESFile)
    file = open(PESFile,"w")
    # header
    file.write("# Potential Energy Surface data\n")
    file.write("# r/A Energy/Hartree \n")
    # write value out
    for iLine in range(Nr):
        file.write(str(rList[iLine])+" "+str(PESList[iLine])+"\n")

    return rList, PESList

# Read single point energy, can only be used for single point energy calculation (N/A for PES calculation)
def ReadSPE(OutputFile):
    # generate SPE file name, in .spe
    SPEFile = OutputFile.split(".out")[0]+"_SPE.spe"
    with open(OutputFile,"r") as f:
            lines = f.read().split("\n")
            for i,line in enumerate(lines):
                if "FINAL SINGLE POINT ENERGY" in line:
                    words =line.split()
                    if len(words)>6:
                        SPE = float(words[6])
                    else:
                        SPE = float(words[4])
    if os.path.exists(SPEFile):
        os.remove(SPEFile)
    file = open(SPEFile,"w")
    # header
    file.write("# Single point energy data\n")
    file.write("# Energy/Hartree \n")
    file.write(str(SPE)+"\n")
    return SPE

# extract PES from PES file
def GetPES(PESFile):
    rList = []
    PESList = []
    with open(PESFile,"r") as f:
        lines = f.read().split("\n")
        for line in lines[2:-1]:
            rList.append(float(line.split(" ")[0]))
            PESList.append(float(line.split(" ")[1]))
    return rList, PESList

def GetSPE(SPEFile):
    with open(SPEFile,"r") as f:
        lines = f.read().split("\n")
        SPE = float(lines[2])
    return SPE

class AtomInfo(object):
 # get a list of atom types and coords
 def __init__(self, mol):
    mol.make3D()
    self.type   = [] 
    self.coord  = []
    self.mass   = []
    self.M      = 0.
    self.NAtoms = len(mol.atoms)
    for iAtom in range(self.NAtoms):
        atom = mol.atoms[iAtom]
        self.type.append(atom.type.split("3")[0])
        x = atom.coords[0]
        y = atom.coords[1]
        z = atom.coords[2]
        self.coord.append([x,y,z])
        self.mass.append(atom.atomicmass)
        self.M += atom.atomicmass
    self.GetCoM()
    self.SetCoMOrigin()

 def GetCoM(self):
    # get center of mass
    moment = [0., 0., 0.]
    for iAtom in range(len(self.type)):
        for jCoord in range(3):
            moment[jCoord] += self.coord[iAtom][jCoord]*self.mass[iAtom]
    self.CoM = [iMoment / self.M for iMoment in moment]

 def SetCoMOrigin(self):
    # set the molecule center of mass to the origin 
    for iAtom in range(len(self.type)):
        for jCoord in range(3):
            self.coord[iAtom][jCoord] -= self.CoM[jCoord]

 def SetCenterOrigin(self):
    # for tetragonal point group Td only, set the center atom to be at the origin
    # using SetToOrigin() based on CoM gives the center coordinate a tiny displacement from origin
    atom_center = self.coord[0].copy()
    for iAtom in range(len(self.type)):
        for jCoord in range(3):
            self.coord[iAtom][jCoord] -= atom_center[jCoord]

 def UnitVector(self, vector):
    # get a unit vector from a given vector
    length = self.Length(vector)
    unitvector = [iVector / length for iVector in vector]
    return unitvector 
 
 def Length(self, vector):
    # get the length of the vector
    length = sqrt( sum(iVector*iVector for iVector in vector) )
    return length

 def Reflection(self, axis):
    # reflection about the plane that is perpandicular to the given axis and contains the origin
    for iAtom in range(len(self.type)):
        if axis == "x":
            self.coord[0] = -self.coord[0]
        if axis == "y":
            self.coord[1] = -self.coord[1]
        if axis == "z":
            self.coord[2] = -self.coord[2]

 def TdOrient(self, config = "A"):
    self.SetCenterOrigin()
    UVec12 = [self.UnitVector(self.coord[1]), self.UnitVector(self.coord[2])]
    if config == "A":
        # the first and initial orientation, A (CH4 as an example) 
        #    H2                  z        
        #     \                  ^       
        #      \                 |       
        #       \                |    
        #       .C------- H1     |----->y
        #     .'/|              /
        #  H3' / |             V
        #     /__|             x
        #      H4             
        r , rmsd = R.align_vectors([[0.,1.,0.],[0., -1./3., sin(acos(-1./3.))]], UVec12)

    elif config == "B":
        # Orientation B (CH4 as an example)
        #      H3                                 
        #     \¯¯|          
        #  H4. \ |           
        #     '.\|          
        #       'C------- H1
        #       /
        #      / 
        #     /
        #    H2
        r , rmsd = R.align_vectors([[0.,1.,0.],[0., -1./3.,-sin(acos(-1./3.))]], UVec12)

    elif config == "C":
        # Orientation C (CH4 as an example)
        #                H2
        #               /
        #              /
        #             /
        #  H1 -------C.
        #            |\'.
        #            | \ 'H4
        #            |__\
        #             H3
        r , rmsd = R.align_vectors([[0.,-1.,0.],[0., 1./3.,sin(acos(-1./3.))]], UVec12)

    elif config == "D":
        # Orientation D (CH4 as an example)
        #             H4
        #            |¯¯/
        #            | / .H3
        #            |/.'
        #  H1 -------C'
        #             \
        #              \ 
        #               \
        #                H2
        r , rmsd = R.align_vectors([[0.,-1.,0.],[0., 1./3.,-sin(acos(-1./3.))]], UVec12)
    
    elif config == "E":
        # Orientation E (CH4 as an example)
        #    H1
        #     \
        #      \       
        #       \ ,,.--'' H3
        #        C________
        #       / ¯`-.__ | H4
        #      /        ¯`
        #     /
        #    H2
        r , rmsd = R.align_vectors([[0., -cos(acos(-1./3.)/2.),  sin(acos(-1./3.)/2.)],\
                                    [0., -cos(acos(-1./3.)/2.), -sin(acos(-1./3.)/2.)]], UVec12)

    elif config == "F":
        # Orientation F (CH4 as an example)
        #                H1
        #               /
        #              /    
        #  H4 ''--... /
        #    ________C
        # H3 | __.-'¯ \
        #    '¯        \ 
        #               \
        #                H2
        r , rmsd = R.align_vectors([[0., cos(acos(-1./3.)/2.), sin(acos(-1./3.)/2.)],\
                                    [0., cos(acos(-1./3.)/2.),-sin(acos(-1./3.)/2.)]], UVec12)
    else:
        raise Exception("config must be A or B or C or D or E or F")

    for iAtom in range(len(self.type)):
        self.coord[iAtom] = r.apply(self.coord[iAtom]) 

class const(object):
    # useful constants
    def __init__(self):
        # ORCA 4.2.1 manual
        self.Hartree2eV = 27.2113834
        self.eV2cm = 8065.54477
        # NIST CODATA 2018
        self.eV2K       = 1.160451812e4
        self.Hartree2K  = self.Hartree2eV*self.eV2K
        # ORCA 4.2.1 manual
        self.Bohr2A     = 0.5291772083
        # Boltzmann in eV
        self.Boltzmann  = 8.617333262e-5
        # eV to J
        self.eV2J       = 1.602176634e-19
        # Avogadro's const
        self.NA         = 6.02214076e23
        # NIST CCCBDB
        self.Hartree2kcal = 627.5

#--------------------------------------------------------------
# Old ORCA_SAFT class

class ORCA_SAFT(object):
 def __init__(self,smiles,method="HFLD_pVDZ.txt",path_orca="orca",path_multiwfn="multiwfn"):
    
    # Create a Molecule object in pybel
    mol = pybel.readstring("smi",smiles)

    # Set it up in 3D coordinates
    mol.make3D()

    self.smiles = smiles
    self.method = method
    self.path_orca = path_orca
    self.path_multiwfn = path_multiwfn
    self.AtomInfo = []
    self.AtomInfo.append(deepcopy(AtomInfo(mol)))
    self.AtomInfo.append(deepcopy(AtomInfo(mol)))
    AI = self.AtomInfo
    #----------------------------------------------------------------------
    # need a handle to input the orientations from job submission 
    TdOriList = ["AA", "AB", "AC", "AD", "CA", "CB", "EE", "EF", "AE", "AF", "FA", "CF"]
    #----------------------------------------------------------------------

    # Making pure input file
    if os.path.exists(smiles+"_pure.inp"):
        os.remove(smiles+"_pure.inp")
    file = open(smiles+"_pure.inp","w")
    with open(method) as f:
        lines = f.read().split("\n")
        for i in range(len(lines)):
            file.write(lines[i]+"\n")
    file.write("\n")
    file.write("* xyz 0 1\n")

    # Obtain the atom coordinates for a molecule
    # N.B. This may need some tweaking as this is not often reproducible
    
    # the old way to assign coordinates:
    # for i in range(len(mol.atoms)):
    #     atom = mol.atoms[i]
    #     atom_type = atom.type.split("3")[0]
    #     x = atom.coords[0]
    #     y = atom.coords[1]
    #     z = atom.coords[2]
    #     file.write("  "+atom_type+"(1)\t"+str(x)+"\t"+str(y)+"\t"+str(z)+"\n")

    for i in range(AI[0].NAtoms):
        file.write("  "+AI[0].type[i]+"(1)\t"+str(AI[0].coord[i][0])+"\t"+str(AI[0].coord[i][1])+"\t"+str(AI[0].coord[i][2])+"\n")
    file.write("*\n\n\n")    
    file.close()
    
    #----------------------------------------------------------------------
    # need modify
    # Add step here to determine Vol and good guess for sigma from multiwfn
    if os.path.exists(smiles+"_opt.inp"):
        os.remove(smiles+"_opt.inp")
    file = open(smiles+"_opt.inp","w")
    with open(method) as f:
        lines = f.read().split("\n")
        for i in range(len(lines)):
            if "!" in lines[i]:
                file.write("! RKS B3LYP D3BJ Opt"+"\n")
    file.write("\n")
    file.write("* xyz 0 1\n")
    for i in range(AI[0].NAtoms):
        file.write("  "+AI[0].type[i]+"(1)\t"+str(AI[0].coord[i][0])+"\t"+str(AI[0].coord[i][1])+"\t"+str(AI[0].coord[i][2])+"\n")
    for i in range(AI[0].NAtoms):
        file.write("  "+AI[0].type[i]+"(1)\t"+str(AI[1].coord[i][0])+"\t"+str(AI[1].coord[i][1])+"\t"+str(AI[1].coord[i][2])+"\n")
    file.write("*\n\n\n")    
    file.close()

    print("orca "+smiles+"_opt.inp>"+smiles+"_opt.out")
    subprocess.run(path_orca+" "+smiles+"_opt.inp>"+smiles+"_opt.out", shell=True)
    
    x = []
    y = []
    z = []
    with open(smiles+"_opt.xyz","r") as f:
        lines = f.read().split("\n")
        for i,line in enumerate(lines[2:len(lines)-1]):
            words =line.split()
            x.append(float(words[1]))
            y.append(float(words[2]))
            z.append(float(words[3]))
    r0 = ((x[0]-x[1])**2+(y[0]-y[1])**2+(z[0]-z[1])**2)**(0.5)
    r = np.linspace(0.95*r0,2*r0,50)
    self.r = r
    #----------------------------------------------------------------------

    # Creating dimer input file
    if os.path.exists(smiles+"_dimer.inp"):
        os.remove(smiles+"_dimer.inp")
    file = open(smiles+"_dimer.inp","w")
    with open(method) as f:
        lines = f.read().split("\n")
        for i in range(len(lines)):
            file.write(lines[i]+"\n")
    file.write("\n")
    file.write("%coords\n")
    file.write("CTyp\txyz\n")
    file.write("Charge\t0\n")
    file.write("Mult\t1\n")
    file.write("pardef\n")
    file.write("\tr = "+str(0.95*r0)+", "+str(2*r0)+", "+str(int(50))+";"+"\n")
    file.write("end\n")
    file.write("coords\n")
    
    # dimer file
    # Molecule 1
    for i in range(AI.NAtoms):
        file.write("  "+AI.type[i]+"(1)\t"+str(AI[0].coord[i][0])+"\t"+str(AI[0].coord[i][1])+"\t"+str(AI[0].coord[i][2])+"\n")
        
    # Molecule 2
    for i in range(AI.NAtoms):
        file.write("  "+AI.type[i]+"(2)\t"+str(AI[1].coord[i][0])+"+{r}"+"\t"+str(AI[1].coord[i][1])+"\t"+str(AI[1].coord[i][2])+"\n")

    file.write("end\n")
    file.write("end\n")
    file.write("%method\n")
    file.write("ScanGuess MORead\n")
    file.write("end\n")
    file.write("\n\n\n")    
    file.close()

 def MolVol(self,cutoff=0.001):
    # Attempt to obtain V_QC from multiwfn
    smiles = self.smiles
    path_multiwfn = self.path_multiwfn

    # Create input file
    if os.path.exists(smiles+"_pure_wfn.txt"):
        os.remove(smiles+"_pure_wfn.txt")
    file = open(smiles+"_pure_wfn.txt","w")
    file.write("12\n")
    file.write("1\n")
    file.write("1\n")
    file.write(str(cutoff)+"\n")
    file.write("0\n")
    file.close()

    # Run input file
    subprocess.run(path_multiwfn+" "+smiles+"_pure_wfn.txt>"+smiles+"_VdW_vol.txt", shell=True)

    # Read output file
    with open(smiles+"_VdW_vol.txt","r") as f:
        lines = f.read().split("\n")
        for j,line in enumerate(lines):
            if "Volume:" in line:
                words = line.split()
                Vol = float(words[4])
    return Vol
 def Potential(self):
    smiles = self.smiles
    r      = self.r
    path_orca = self.path_orca
    # Run pure input file
    print("orca "+smiles+"_pure.inp>"+smiles+"_pure.out")
    subprocess.run(path_orca+" "+smiles+"_pure.inp>"+smiles+"_pure.out", shell=True)
    
    with open(smiles+"_pure.out","r") as f:
        lines = f.read().split("\n")
        for i,line in enumerate(lines):
            if "FINAL SINGLE POINT ENERGY" in line:
                words =line.split()
                if len(words)>6:
                    V_pure = float(words[6])
                else:
                    V_pure = float(words[4])
    
    # Run dimer input file
    print("orca "+smiles+"_dimer.inp>"+smiles+"_dimer.out")
    subprocess.run(path_orca+" "+smiles+"_dimer.inp>"+smiles+"_dimer.out", shell=True)
    
    V = np.zeros(np.shape(r))
    with open(smiles+"_dimer.out","r") as f:
        lines = f.read().split("\n")
        i=0
        for j,line in enumerate(lines):
            if "FINAL SINGLE POINT ENERGY" in line:
                words =line.split()
                if len(words)>6:
                    V[i] = (float(words[6])-2*V_pure)*315775.082
                else:
                    V[i] = (float(words[4])-2*V_pure)*315775.082
                i+=1
    return V
 def Clean(self):
    smiles = self.smiles

    # Remove files that aren't necessary
    for filename in glob.glob("./*"):
        if smiles in filename:
            # Don't remove input and output files
            if ("out" not in filename) and ("inp" not in filename):
                os.remove(filename)
 def SAFT_Params(self):
    
    # Obtain the QC potential
    V = self.Potential()
    r = self.r
    self.V = V

    # Obtain sigma and epsilon from QC
    epsilon = -np.amin(V)
    sigma = interp1d(V[0:np.where(V==-epsilon)[0][0]],r[0:np.where(V==-epsilon)[0][0]])(0)
    print("epsilon  = "+str(epsilon)+" K")
    print("sigma    = "+str(sigma)+" angstrom")

    self.epsilon = epsilon
    self.sigma = sigma

    # Optimise for lambda_R
    cons = ({'type': 'ineq', 'fun': lambda x:  x[0]-6})
    res = minimize(self.Obj,[12],args=(r/sigma,V/epsilon),constraints=cons)
    lambda_R = res.x[0]
    print("lambda_A = 6")
    print("lambda_R = "+str(lambda_R))

    self.lambda_A = 6
    self.lambda_R = lambda_R
    return [epsilon,sigma,lambda_R]

 def Obj(self,x,*arg):
    # Objective function for fitting lambda_R
    r,V = arg
    V_Mie = (x[0]/(x[0]-6))*(x[0]/6)**(6/(x[0]-6))*((1/r)**x[0]-(1/r)**6)
    return np.sum((((V-V_Mie)/V)**2)*np.exp(-1-V))

 def Plotting(self):
    epsilon = self.epsilon
    sigma   = self.sigma
    lambda_A = self.lambda_A
    lambda_R = self.lambda_R
    
    r     = self.r
    V     = self.V

    # Obtain Mie potential
    V_Mie = epsilon*(lambda_R/(lambda_R-lambda_A))*(lambda_R/lambda_A)**(lambda_A/(lambda_R-lambda_A))*((sigma/r)**lambda_R-(sigma/r)**lambda_A)

    plt.rc('font',family='serif',size=16)
    plt.rc('text',usetex=True)
    plt.rc('xtick',labelsize='small')
    plt.rc('ytick',labelsize='small')
    fig = plt.figure(figsize=(7,5.25))
    ax = fig.add_subplot(1,1,1)
    ax.plot(r,V,color='r',ls='solid',linewidth=2,label=r'QC')
    ax.plot(r,V_Mie,color='b',ls='solid',linewidth=2,label=r'Mie')
    ax.set_xlabel(r'$r / \AA$')
    ax.set_ylabel(r'$\phi^\mathrm{Mie} /$ K')
    ax.legend(frameon=False)
    # ax.set_ylim(-1.1, 2)
    # ax.set_xlim(0.6,1.9)
    plt.show()

# old ORCA_SAFT class end
# ---------------------------------------------------------------

