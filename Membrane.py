import numpy as np
from Bio.PDB.PDBParser import PDBParser

class MembranePlacer(object):
    def __init__(self,helices,structure,normal):
        self.coords={}
        self.normal=normal
        self.uppermembrane=np.array([0,0,0])
        self.lowermembrane=np.array([0,0,0])
        self.helices=helices
        self.structure=structure

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
        return averagemiddle-10, averagemiddle+10


    #def shiftMembrane(self):

    def score(self):
        hydrophobic_residues = ["Phe", "Gly", "Ile", "Leu", "Mis", "Val", "Trp", "Thr"]
        hydrophilic_residues = ["Ala", "Cys", "Asp", "Glu", "His", "Lys", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr"]
        hphobcount=0
        hphilcount=0
        for residue in self.structure.get_residues():
            try:
                atom = residue['CA']
                coordinates = atom.get_coord()
                up=coordinates-self.uppermembrane
                down=coordinates-self.lowermembrane
                first=up.dot(self.uppermembrane)
                second=down.dot(self.lowermembrane)
                if (first<0 and second>0) or (first>0 and second<0):
                    if hydrophilic_residues.__contains__(residue):
                        hphilcount+=1
                    else:
                        hphobcount+=1
            except KeyError:
                print("no CA")
        print(hphilcount)
        print(hphobcount)


    def placeMembrane(self):
        self.lowermembrane= self.firstGuess()[0]
        self.uppermembrane= self.firstGuess()[1]
        self.score()
        return self.firstGuess()