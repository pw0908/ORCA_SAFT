import os
import numpy as np
from math import exp, sqrt, acos, cos, sin
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation as R
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from openbabel import pybel
from copy import deepcopy
import const
from shutil import copy

class GenerateInput(object):
 """
 Generate input files for ORCA calculations. 
 
 These input files are for the calculations of single molecule and dimer with different orientations.
 
 Parameters
 ----------
 smiles: string
     a SMILES string for molecule
 method: string
     name of header files for generating ORCA input files (without `.txt`); these input header files 
     includes all information about basis sets and quantum chemistry methods, example files can be 
     found in `ORCA_header` folder. 
     
 Returns
 -------
 no return
 """
 def __init__(self,smiles,method="HFLD_pVDZ"):
    
    # Create a Molecule object in pybel
    mol = pybel.readstring("smi",smiles)

    # Set it up in 3D coordinates
    mol.make3D()

    self.smiles = smiles
    self.method = method
    self.AtomInfo = []
    self.AtomInfo.append(deepcopy(AtomInfo(mol)))
 
 def InputPure(self, mode="single", suffix="pure", nprocs=8):
    """
    Generate calculation input for a single molecule, without counterpoise correction
    
    Parameters
    ----------
    mode: string, "single" or "parallel"
        mode of calculations, 
        if "single", the generated file will be the input for a single processor serial calculation
        if "parallel", the generated file will be the input for a parallel calculation 
        if "parallel" mode is chosen, the ncpus variable need to be specified (default = 8), which is
        the number of cpus for calculation
    suffix: string
        suffix for the input file name, the input file name has the following structure:
        <SMILES of the species>_<suffix>_<method>.inp
    nprocs: int > 0
        number of processors to be used in parallel computations.
        
    Returns
    -------
    no return, generate an input file for single molecule calculation without counterpoise correction
    """
    AI = self.AtomInfo
    # generate input file name string
    PureFileName = self.smiles+"_"+suffix+"_"+self.method+".inp"
    # test if file exists, if so, remove it and generate a new file.
    if os.path.exists(PureFileName):
        os.remove(PureFileName)
    file = open(PureFileName,"w")

    # copy the header file info
    with open("./ORCA_header/"+self.method+".txt") as f:
        lines = f.read().split("\n")
        for i in range(len(lines)):
            file.write(lines[i]+"\n")
    
    # specify processor number
    if mode == "parallel":
        file.write("%pal nprocs "+str(int(nprocs))+" end\n")

    file.write("\n")
    
    # write coordinates from AtomInfo
    file.write("* xyz 0 1\n")
    for iAtom in range(AI[0].NAtoms):
        file.write("  "+AI[0].type[iAtom]+"(1)\t"+str(AI[0].coord[iAtom][0])+
        "\t"+str(AI[0].coord[iAtom][1])+"\t"+str(AI[0].coord[iAtom][2])+"\n")
    file.write("*\n\n\n")    
    file.close()
 
 def InputPureCP(self, rval, molecule=1, mode = "single", suffix = "pure", 
    nprocs = 8, group = None, ori = None):
    """
    Generate calculation input for a single molecule, with counterpoise correction.
    
    By now it supports calculations for a single atom, or molecules with tetragonal (Td) symmetry.
    
    Parameters
    ----------
    rval: float > 0, or a list of floats > 0
        the center of mass distance of two species
        can be either a float number or a list of float numbers
    molecule: int, 1 or 2
        specify energy of which species is calculated, 1 corresponds to the species with its center of
        mass located at the origin
    mode: string, "single" or "parallel"
        mode of calculations, 
        if "single", the generated file will be the input for a single processor serial calculation
        if "parallel", the generated file will be the input for a parallel calculation 
        if "parallel" mode is chosen, the ncpus variable need to be specified (default = 8), which is
        the number of cpus for calculation
    suffix: string
        suffix for the input file name, the input file name has the following structure:
        for atom: 
        <SMILES of the species>_<suffix>_CP_<method>.inp
        for molecule with Td symmetry: 
        <SMILES of the species>_<suffix>_CP_<method>_<group>_<ori>.inp
    nprocs: int > 0
        number of processors to be used in parallel computations.
    group: None or string "Td"
        The symmetry group for species; if the calculation is performed for atom, then the default None
        is used, if the molecule has Td symmetry (which is the only supported group), then group = "Td".
    ori: string in {"AA", "AB", "AC", "AD", "CA", "CB", "EE", "EF", "AE", "AF", "FA", "CF"}
        Only supports for Td symmetry, specify orientations of the species and the ghost species. The 
        first letter is for the species 1 and the second letter is for the species 2. See readme for 
        the definition of these orientations.
    Returns
    -------
    no return, generate an input file for single molecule calculation without counterpoise correction
    """
    if molecule not in {1,2}:
        raise Exception("molecule need to be either 1 or 2")
    # extract species information
    AI = self.AtomInfo
    # so far considered only Td geometry 
    if group == "Td":
        AI[0].TdOrient(ori[0])
        AI[1].TdOrient(ori[1])
    # possibly write a function in future to systematically check the two-body
    # orientations in a systematic way.
    else:
        pass
    # create file name and initialize file
    DimerFileName = self.smiles+"_"+suffix+"_CP_"+str(molecule)+"_"+self.method
    if None not in {group, ori}:
        DimerFileName += ("_"+group+"_"+ori)
    DimerFileName += ".inp"

    # copy the header file
    if os.path.exists(DimerFileName):
        os.remove(DimerFileName)
    file = open(DimerFileName,"w")
    # preset methods
    with open("./ORCA_header/"+self.method+".txt") as f:
        lines = f.read().split("\n")
        for i in range(len(lines)):
            file.write(lines[i]+"\n")

    if mode == "parallel":
        file.write("%pal nprocs "+str(int(nprocs))+" end\n")
    
    # write coordinate block
    file.write("\n")
    file.write("%coords\n")
    file.write("CTyp\txyz\n")
    file.write("Charge\t0\n")
    file.write("Mult\t1\n")

    # define a parameter for varying r
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
    
    # dimer position, counterpoise correction will make one of the molecule to 
    # be the "ghost" - with only the basis sets left there, the electrons and
    # ions are not.
    # for the case where real molecule is the molecule 1
    if molecule == 1:
        # Molecule 1
        for iAtom in range(AI[0].NAtoms):
            file.write("  "+AI[0].type[iAtom]+"\t"\
            +str(AI[0].coord[iAtom][0])+"\t"\
            +str(AI[0].coord[iAtom][1])+"\t"\
            +str(AI[0].coord[iAtom][2])+"\n") 
        # Molecule 2 (ghost)
        for iAtom in range(AI[1].NAtoms):
            file.write("  "+AI[1].type[iAtom]+":\t"\
            +str(AI[1].coord[iAtom][0])+"\t"\
            +"{"+ str(AI[1].coord[iAtom][1]) +"+r}\t"\
            +str(AI[1].coord[iAtom][2])+"\n")
    # for the case where real molecule is the molecule 2
    elif molecule == 2:
        # Molecule 1 (ghost)
        for iAtom in range(AI[0].NAtoms):
            file.write("  "+AI[0].type[iAtom]+":\t"\
            +str(AI[0].coord[iAtom][0])+"\t"\
            +str(AI[0].coord[iAtom][1])+"\t"\
            +str(AI[0].coord[iAtom][2])+"\n")
        # Molecule 2 
        for iAtom in range(AI[1].NAtoms):
            file.write("  "+AI[1].type[iAtom]+"\t"\
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

 def InputDimer(self, rval, mode = "single", suffix = "dimer", nprocs = 8,\
    group = None, ori = None):
    """
    Generate calculation input for a pair of species.
   
    By now it supports calculations for a single atom, or molecules with tetragonal (Td) symmetry.
    
    Parameters
    ----------
    rval: float > 0, or a list of floats > 0
        the center of mass distance of two species
        can be either a float number or a list of float numbers
    mode: string, "single" or "parallel"
        mode of calculations, 
        if "single", the generated file will be the input for a single processor serial calculation
        if "parallel", the generated file will be the input for a parallel calculation 
        if "parallel" mode is chosen, the ncpus variable need to be specified (default = 8), which is
        the number of cpus for calculation
    suffix: string
        suffix for the input file name, the input file name has the following structure:
        for atoms:
        <SMILES of the species>_<suffix>_<method>.inp
        for molecules:
        <SMILES of the species>_<suffix>_<method>_<group>_<ori>.inp
    nprocs: int > 0
        number of processors to be used in parallel computations.
    group: None or string "Td"
        The symmetry group for species; if the calculation is performed for atom, then the default None
        is used, if the molecule has Td symmetry (which is the only supported group), then group = "Td".
    ori: string in {"AA", "AB", "AC", "AD", "CA", "CB", "EE", "EF", "AE", "AF", "FA", "CF"}
        Only supports for Td symmetry, specify orientations of the species. The first letter is for the
        species 1 and the second letter is for species 2. See readme for the definition of these 
        orientations.
    Returns
    -------
    no return, generate an input file for a pair of species calculation
    """
    if orcaproc is None:
        nprocs = ncpus
    else:
        nprocs = orcaproc
    # extract species information
    AI = self.AtomInfo
    # so far considered only Td geometry 
    if group == "Td":
        AI[0].TdOrient(ori[0])
        AI[1].TdOrient(ori[1])
    # possibly write a function in future to systematically check the two-
    # body orientations in a systematic way.
    else:
        pass
    # generate dimer file name
    DimerFileName = self.smiles+"_"+suffix+"_"+self.method
    if None not in {group, ori}:
        DimerFileName += ("_"+group+"_"+ori)
    DimerFileName += ".inp"
    if os.path.exists(DimerFileName):
        os.remove(DimerFileName)
    file = open(DimerFileName,"w")
    # preset methods
    with open("./ORCA_header/"+self.method+".txt") as f:
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

def ReadPES(*, Nr, OutputFile):
    """
    Read the potential energy surface (PES) data from the output file as specified in OutputFile,
    also generate a file which is a list of r values and energies.
    
    Suitable for outputs of single-species calculations with counterpoise correction, and pair 
    calculations.
    
    A note on energy types: in some calculations, there are different types of energies, for HFLD, there 
    are "'Actual Energy'", "SCF Energy" and "MDCI Energy". The "'Actual Energy'" is the same as 
    "SCF Energy", but they are not the energy we are looking for. We need "MDCI Energy" because the CI 
    part between the two fragments are added to the SCF energy. 
    
    Parameters
    ----------
    Nr: int > 0
        the # of r calculated
    OutputFile: string
        name of the output files, with the suffix ".out"
        
    Returns
    -------
    rList, PESList: a list of r values and energy values 
    also generate a file of these values, with file name
    <OUtputFile without .out>_PES.pes
    """
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
            # find the correct line for PES data and read it
            if "FINAL SINGLE POINT ENERGY" in line:
                words = line.split()
                if len(words)>6:
                    PESList.append(float(words[6]))
                else:
                    PESList.append(float(words[4]))
                iPES += 1
            # find the r value and read it
            if "R  :" in line:
                words = line.split("R  :")
                rList.append(float(words[1]))
                ir += 1
    if (iPES != Nr):
        raise Exception("Either the # of r or the # of energy does not match \
            the desired #, please check outputs")

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

def ReadSPE(OutputFile):
    """
    Read single point energy, can only be used for single point energy calculations
    
    Parameters
    ----------
    OutputFile: string
        name of the output files, with the suffix ".out"
        
    Returns
    -------
    SPE: a single point energy
    """
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
    # write the SPE in a file
    file.write("# Single point energy data\n")
    file.write("# Energy/Hartree \n")
    file.write(str(SPE)+"\n")
    return SPE

def MergePES(PESin1, PESin2, PESNew = None):
    """
    merge two PES files
    
    Parameters
    ----------
    PESin1: string
        name of the 1st PES file, with the file extension name
    PESin2: string
        name of the 2nd PES file, with the file extension name
    PESNew: string or None
        name of the new PES file, with the file extension name 
        if it is None, the 2nd file will be added onto the 1st file
        
    Returns
    -------
    no return, either the 1st file is appended with the 2nd file, or a new file is generated
    """
    if os.path.isfile('./'+PESin1) + os.path.isfile('./'+PESin2) != 2:
        raise Exception("input file not found")
    elif PESNew is not None:
        if (os.path.isfile('./'+PESNew)) is True:
            raise Exception(PESNew+" already exists")
    
    with open(PESin2,"r") as f2:
        lines2 = f2.read().split("\n")[2:]
    
    if PESNew is None:
        f1w = open(PESin1,"a")
    else:
        with open(PESin1, "r") as f1:
            lines1 = f1.read()
        f1w = open(PESNew,"w")
        for line1 in lines1:
            f1w.write(line1)
    
    for line2 in lines2:
        f1w.write(line2+"\n")

    f1w.close()

def GetPES(PESFile):
    """
    extract PES from PES file
    
    Parameters
    ----------
    PESFile: string
        PES file name
        
    Returns
    -------
    rList, PESList: a list of r values and energy values 
    """
    rList = []
    PESList = []
    with open(PESFile,"r") as f:
        lines = f.read().split("\n")
        lines = [jline for jline in lines if jline != '']
        for line in lines[2:-1]:
            rList.append(float(line.split(" ")[0]))
            PESList.append(float(line.split(" ")[1]))
    # remove duplicated r and PES values
    RemoveDup(rList)
    RemoveDup(PESList)
    if len(rList) != len(PESList):
        raise Exception("Inconsistent data in PES file, please check")
    return rList, PESList

def RemoveDup(seq):
    """
    remove duplicated entries in a list, if any
    
    Parameters
    ----------
    seq: list
        a list with possible duplicates
    
    Returns
    -------
    a list with no duplicate
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def GetSPE(SPEFile):
    """
    a function to extract SPE from the file
    
    Parameters
    ----------
    SPEFile: string
        the file name for the SPE file
    
    Returns
    -------
    SPE value
    """
    with open(SPEFile,"r") as f:
        lines = f.read().split("\n")
        SPE = float(lines[2])
    return SPE

class AtomInfo(object):
 """
 get a list of atom types and coordinations, with a collection of helper functions
 
 Parameters
 ----------
 mol: Pybel molecule object
     molecule object for the species of interests
 
 Returns
 -------
 no return
 """
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
    """
    obtain the center of mass of the species
    """
    # get center of mass
    moment = [0., 0., 0.]
    for iAtom in range(len(self.type)):
        for jCoord in range(3):
            moment[jCoord] += self.coord[iAtom][jCoord]*self.mass[iAtom]
    self.CoM = [iMoment / self.M for iMoment in moment]

 def SetCoMOrigin(self):
    """
    set the molecule center of mass to the origin 
    """
    for iAtom in range(len(self.type)):
        for jCoord in range(3):
            self.coord[iAtom][jCoord] -= self.CoM[jCoord]

 def SetCenterOrigin(self):
    """
    for tetragonal point group Td only, set the center atom to be at the origin
    using SetToOrigin() based on CoM gives the center coordinate a tiny displacement from origin
    """
    atom_center = self.coord[0].copy()
    for iAtom in range(len(self.type)):
        for jCoord in range(3):
            self.coord[iAtom][jCoord] -= atom_center[jCoord]

 def UnitVector(self, vector):
    """
    get a unit vector from a given vector
    """
    length = self.Length(vector)
    unitvector = [iVector / length for iVector in vector]
    return unitvector 
 
 def Length(self, vector):
    """
    get the length of a vector
    """
    length = sqrt( sum(iVector*iVector for iVector in vector) )
    return length

 def Reflection(self, axis):
    """
    reflection about the plane that is perpandicular to the given axis and contains the origin
    """
    for iAtom in range(len(self.type)):
        if axis == "x":
            self.coord[0] = -self.coord[0]
        if axis == "y":
            self.coord[1] = -self.coord[1]
        if axis == "z":
            self.coord[2] = -self.coord[2]

 def TdOrient(self, config = "A"):
    """
    for a molecule with Td symmetry, orient the molecule in a given orientation
   
    Parameters
    ----------
    config: a character {"A", "B", "C", "D", "E", "F"}
        possible orientation configurations; for the meaning of these, see readme
    """
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

class Potential(object):
 """
 A class that takes numerical r and V values as inputs and obtain all the Mie parameters.
 
 Parameters
 ----------
 r: list of float > 0
     the inter species distance
 V: list of float
     a list of potential values, in Hartree
 
 Returns
 -------
 no return
 """
 def __init__(self, r, V):
    self.r = r 
    # now energy is in K
    self.V = V * const.Hartree2K

 def SAFT_Params(self):
    """
    Obtain SAFT parameters for potentials obtained by QC calculations with no analytic form
    
    Parameters
    ----------
    no parameters
    
    Returns
    -------
    a list of epsilon, sigma, and lambda_R values
    also the following attributes are generated: epsilon, sigma, lambda_A, lambda_R
    """
    r = self.r 
    V = self.V 
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
