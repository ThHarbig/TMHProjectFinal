__author__ = "Kevin Menden"
__date__ = '26.04.2016'

import sys, os, glob
import operator
import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser

class PDBHelixParser(object):

    def __init__(self, file):
        self.file = file
        self.totalLength = 0
        self.alphaHelix = 0
        self.helices = []
        self.structure = None
        self.proteinHelixSequences = []
        self.alphaHelixPositions = []
        self.aminoAcidDict = {"ALA": "A", "GLY": "G", "PHE": "F", "ILE": "I", "MET": "M", "LEU": "L", "PRO": "P", "VAL": "V",
                          "ASP": "D", "GLU": "E", "LYS": "K", "ARG": "R", "SER": "S", "THR": "T", "TYR": "Y", "HIS": "H",
                          "CYS": "C", "ASN": "N", "GLN": "Q", "TRP": "W", "MSE": "M" }

        self.sortedHelix = []
        self.backboneDistances = []
        self.backboneDistHelix = []


    def getResiduesFromChain(self, chain, start, end):
        """
        Extract the amino acids given a chain and the end and start positions
        :param chain: the chain
        :param start: start of the helix
        :param end: end of the helix
        :return: array containing all residues of the helix
        """
        helix_residues = []
        for residue in chain.get_residues():
            if (residue.id[1] >= start and residue.id[1] <= end and residue.get_resname() in self.aminoAcidDict.keys()):
                helix_residues.append(residue)

        return helix_residues

    def getResiduesFromList(self, residues, start, end):
        """
        :param residues:
        :param start:
        :param end:
        :return:
        """
        helix_residues = []
        for res in residues:
            if (res.id[1] >= start and res.id[1] <= end and res.get_resname() in self.aminoAcidDict.keys()):
                helix_residues.append(res)

        return helix_residues



    def parsePDBInformation(self, file):
        """
        Parses a single pdb file and counts the residues in
        helices, sheets and the total length of the protein
        """
        helices = []
        helixSequences = []
        f = open(self.file)
        line = f.readline()
        while line:
            #If HELIX, check for type and add length and positions
            if line.startswith("HELIX"):
                start = int(line[21:25].replace(" " ,""))
                end = int(line[33:37].replace(" ", ""))
                type = int(line[39:40])
                chain = line[19:20].replace(" ", "")
                currentHelix = (start, end, chain)
                if type == 1:
                    helices.append(currentHelix)
                line = f.readline()
            else:
                line = f.readline()
        f.close()

        # Parse the structure with a PDBParser object
        pdbParser = PDBParser()
        structure = pdbParser.get_structure("pdbFile", self.file)
        self.structure = structure

        # For every helix tuple, extract the residues and store them in helixSequences
        for helix in helices:
            if helix[2] == "":
                residues = structure.get_residues()
                helixSequences.append(self.getResiduesFromList(residues, helix[0], helix[1]))
            chains = structure.get_chains()
            for chain in chains:
                if (chain.get_id() == helix[2]):
                    helixSequences.append(self.getResiduesFromChain(chain, helix[0], helix[1]))
        return helixSequences



    def parse_pdb_file(self):
        """
        Parse a pdb file and return the helix sequences
        """
        if self.file.endswith(".pdb"):
            print("Parsing file: " + str(self.file))
            self.proteinHelixSequences = self.parsePDBInformation(self.file)




