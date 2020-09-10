import subprocess,os,glob
import numpy as np
from math import exp
from scipy.interpolate import interp1d
from scipy.optimize import fsolve,minimize
import matplotlib.pyplot as plt
from openbabel import pybel

class ORCA_SAFT(object):
 def __init__(self,smiles,method="HFLD.txt"):
    mol = pybel.readstring("smi",smiles)
    mol.make3D()
    
    self.smiles = smiles
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
    for i in range(len(mol.atoms)):
        atom = mol.atoms[i]
        atom_type = atom.type.split("3")[0]
        x = atom.coords[0]
        y = atom.coords[1]
        z = atom.coords[2]
        file.write("  "+atom_type+"(1)\t"+str(x)+"\t"+str(y)+"\t"+str(z)+"\n")
    file.write("*\n\n\n")    
    file.close()

    # Add step here to determine Vol and good guess for sigma
    r = np.linspace(3,7,2)
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
            file.write("  "+atom_type+"(2)\t"+str(x)+"\t"+str(y)+"\t"+str(z+r[j])+"\n")
        file.write("*\n\n\n")    
    file.close()
 def MolVol(self,cutoff=0.001):
    smiles = self.smiles
    if os.path.exists(smiles+"_pure_wfn.txt"):
        os.remove(smiles+"_pure_wfn.txt")
    file = open(smiles+"_pure_wfn.txt","w")
    file.write("12\n")
    file.write("1\n")
    file.write("1\n")
    file.write(str(cutoff)+"\n")
    file.write("0\n")
    file.close()
    subprocess.run(["Multiwfn",smiles+"_pure_wfn.txt>"+smiles+"_VdW_vol.txt"], shell=True)

    with open(smiles+"_VdW_vol.txt","r") as f:
        lines = f.read().split("\n")
        for j,line in enumerate(lines):
            words = line.split()
            Vol = float(words[4])
    return Vol
 def Potential(self):
    smiles = self.smiles
    r      = self.r
    # Run pure input file
    print("orca "+smiles+"_pure.inp>"+smiles+"_pure.out")
    subprocess.run(["orca",smiles+"_pure.inp>"+smiles+"_pure.out"], shell=True)
    
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
    subprocess.run(["orca",smiles+"_dimer.inp>"+smiles+"_dimer.out"], shell=True)
    
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
    for filename in glob.glob("./"+smiles+"_pure*"):
        if ("inp" or "out") not in filename:
            os.remove(filename) 
    for filename in glob.glob("./"+smiles+"_dimer*"):
        if ("inp" or "out") not in filename:
            os.remove(filename) 