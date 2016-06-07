__author__ = "Kevin Menden"
__date__ = '01.06.2016'

import numpy as np

class Helix(object):

    def __init__(self, helix_residues):
        self.aminoAcidDict = {"ALA": "A", "GLY": "G", "PHE": "F", "ILE": "I", "MET": "M", "LEU": "L", "PRO": "P", "VAL": "V",
                         "ASP": "D", "GLU": "E", "LYS": "K", "ARG": "R", "SER": "S", "THR": "T", "TYR": "Y", "HIS": "H",
                         "CYS": "C", "ASN": "N", "GLN": "Q", "TRP": "W", "MSE": "M"}
        self.start_point = helix_residues[0]['N'].get_coord()
        self.end_point = helix_residues[-1]['N'].get_coord()
        self.length = helix_residues[-1]['N'] - helix_residues[0]['N']
        self.vector = np.array([self.end_point[0] - self.start_point[0],
                       self.end_point[1] - self.start_point[1],
                       self.end_point[2] - self.start_point[2]])
        self.neg_vector = -1 * self.vector
        self.sequence = ""
        for res in helix_residues:
            self.sequence += self.aminoAcidDict[res.get_resname()]
        self.positions=[]
        for res in helix_residues:
            self.positions.append(res.id[1])