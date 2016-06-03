import numpy as np
from Bio.PDB.PDBParser import PDBParser

class MembranePlacer(object):
    def __init__(self,helices,structure,normal):
        self.coords={}
        self.normal=normal
        self.uppermembrane=np.array()
        self.lowermembrane=np.array()
        self.helices=helices

    def firstGuess(self):
        averagemiddle=np.array([0.,0.,0.])
        for point in self.helices:
            start=point.start_point
            end=point.end_point
            if np.sqrt((start[0]-end[0])**2+(start[1]-end[1])**2+(start[2]-end[2])**2)<50:
                averagemiddle[0]+=((start[0]+end[0])/2)
                averagemiddle[1]+=((start[1]+end[1])/2)
                averagemiddle[2]+=((start[2] + end[2]) / 2)
        for i in range(0,len(averagemiddle)):
            averagemiddle[i]/=len(self.helices)
        return averagemiddle


    #def shiftMembrane(self):

    def score(self):
        for residue in self.residues:
            atom = residue['CA']
            coordinates = atom.get_coord()


    def placeMembrane(self):
        self.firstGuess()
        self.score()