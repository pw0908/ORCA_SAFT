import subprocess,os,glob
import numpy as np
from math import exp, sqrt, acos, cos, sin
from scipy.interpolate import interp1d
from scipy.optimize import fsolve,minimize
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
from openbabel import pybel

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
    for i in range(len(mol.atoms)):
        atom = mol.atoms[i]
        atom_type = atom.type.split("3")[0]
        x = atom.coords[0]
        y = atom.coords[1]
        z = atom.coords[2]
        file.write("  "+atom_type+"(1)\t"+str(x)+"\t"+str(y)+"\t"+str(z)+"\n")
    file.write("*\n\n\n")    
    file.close()

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
    for i in range(len(mol.atoms)):
        atom = mol.atoms[i]
        atom_type = atom.type.split("3")[0]
        x = atom.coords[0]
        y = atom.coords[1]
        z = atom.coords[2]
        file.write("  "+atom_type+"(1)\t"+str(x)+"\t"+str(y)+"\t"+str(z)+"\n")
    for i in range(len(mol.atoms)):
        atom = mol.atoms[i]
        atom_type = atom.type.split("3")[0]
        x = atom.coords[0]
        y = atom.coords[1]
        z = atom.coords[2]
        file.write("  "+atom_type+"(2)\t"+str(x)+"\t"+str(y)+"\t"+str(z+5)+"\n")
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

    # Creating dimer input file
    if os.path.exists(smiles+"_dimer.inp"):
        os.remove(smiles+"_dimer.inp")
    file = open(smiles+"_dimer.inp","w")
    for j in range(len(r)):
        if j!=0:
            file.write("$new_job\n")
        with open(method) as f:
            lines = f.read().split("\n")
            for i in range(len(lines)):
                file.write(lines[i]+"\n")
        file.write("\n")
        file.write("* xyz 0 1\n")

        # Molecule 1
        for i in range(len(mol.atoms)):
            atom = mol.atoms[i]
            atom_type = atom.type.split("3")[0]
            x = atom.coords[0]
            y = atom.coords[1]
            z = atom.coords[2]
            file.write("  "+atom_type+"(1)\t"+str(x)+"\t"+str(y)+"\t"+str(z)+"\n")
        
        # Molecule 2
        for i in range(len(mol.atoms)):
            atom = mol.atoms[i]
            atom_type = atom.type.split("3")[0]
            x = atom.coords[0]
            y = atom.coords[1]
            z = atom.coords[2]
            file.write("  "+atom_type+"(2)\t"+str(x)+"\t"+str(y)+"\t"+str(z+r[j])+"\n")
        file.write("*\n\n\n")    
    file.close()

 def GetAtomInfo(self, mol):
    # get a list of atom types and coords
    atom_type = []
    atom_coord = []
    atom_mass = []
    for iAtom in range(len(mol.atoms)):
        atom = mol.atoms[iAtom]
        atom_type.append(atom.type.split("3")[0])
        x = atom.coords[0]
        y = atom.coords[1]
        z = atom.coords[2]
        atom_coord.append([x,y,z])
        atom_mass.append(atom.atomicmass)
    return atom_type, atom_coord, atom_mass

 def CoM(self, atom_coord, atom_mass):
    # get center of mass
    M = 0.
    moment = [0., 0., 0.]
    for iAtom in range(len(atom_coord)):
        M += atom_mass[iAtom]
        for jCoord in range(3):
            moment[jCoord] += atom_coord[iAtom][jCoord]*atom_mass[iAtom]
    CoM = [iMoment / M for iMoment in moment]
    return CoM
 
 def SetCoMOrigin(self, atom_coord, atom_mass):
    # set the molecule center of mass to the origin 
    CoM = self.CoM(atom_coord, atom_mass)
    for iAtom in range(len(atom_coord)):
        for jCoord in range(3):
            atom_coord[iAtom][jCoord] -= CoM[jCoord]
    return atom_coord

 def SetCenterOrigin(self, atom_coord):
    # for tetragonal point group Td only, set the center atom to be at the origin
    # using SetToOrigin() based on CoM gives the center coordinate a tiny displacement from origin
    atom_center = atom_coord[0].copy()
    for iAtom in range(len(atom_coord)):
        for jCoord in range(3):
            atom_coord[iAtom][jCoord] -= atom_center[jCoord]
    return atom_coord

 def UnitVector(self, vector):
    # get a unit vector from a given vector
    length = Length(vector)
    unitvector = [iVector / length for iVector in vector]
    return unitvector 
 
 def Length(self, vector):
    # get the length of the vector
    length = sqrt( sum(iVector*iVector for iVector in vector) )
    return length

 def Reflection(self, axis, coord):
    # reflection about the plane that is perpandicular to the given axis and contains the origin
    if axis == "x":
        coord[0] = -coord[0]
    if axis == "y":
        coord[1] = -coord[1]
    if axis == "z":
        coord[2] = -coord[2]
    return coord

 def InitializeTdGeom(self, atom_coord):
    # the first geometry
    R.align_vectors([[1.,0.,0.],[-1./3., sin(acos(-1./3.)) , 0.]], [UnitVector(atom_coord[1]), UnitVector(atom_coord[1])])




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
