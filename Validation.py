import os
import sys

from Helix import Helix
from ModuleMethods import is_transmembrane_helix
from PDBHelixParser import PDBHelixParser

__date__ = '08.06.2016'
"""
Validation of the model
"""


def isTransmembraneProtein(file):
    # parser = argparse.ArgumentParser(description="Membrane Plane Finder")
    # parser.add_argument('pdb')
    pdbParser = PDBHelixParser(file)
    pdbParser.parse_pdb_file()
    structure = pdbParser.structure         # The whole structure of the PDB file
    raw_helices = pdbParser.proteinHelixSequences

    # Convert raw helices into Helix objects
    helix_set = []
    for h in raw_helices:
        if len(h) > 0:
            helix_set.append(Helix(h))
    print("Found " + str(len(helix_set)))
    # Predict transmembrane helices
    tmh_set = []
    for h in helix_set:
        if is_transmembrane_helix(h) == 'tm':
            tmh_set.append(h)
    print(len(tmh_set))
    #print(len(tmh_set) / len(helix_set))
    if len(tmh_set) < 3 or len(tmh_set)/len(helix_set) < 0.1:
        return False
    else:
        return True


if __name__ == "__main__":
    tp = 0
    tn = 0
    fp = 0
    fn = 0

    tm = 0
    glob = 0
    test = []
    directory = "validation_set_glob2"
    for subdir, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".pdb"):
                f = directory + "/" + file
                if isTransmembraneProtein(f):
                    tm += 1
                    fp += 1
                else:
                    tn += 1
                    glob += 1

    directory = "validation_set_tm"
    for subdir, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".pdb"):
                f = directory + "/" + file
                if isTransmembraneProtein(f):
                    tm += 1
                    tp += 1
                else:
                    fn += 1
                    glob += 1

    print("Tramsmembrane: " + str(tm))
    print("Globular: " + str(glob))
    sensitvity = tp /(tp + fn)
    print("Sensitivity: " + str(sensitvity))
    specificity = tn / (tn + fp)
    print("Specificity: " + str(specificity))
    accuray = (tp + tn) / (tp + tn + fp + fn)
    print("Accuracy: " + str(accuray))
    print("TP :" + str(tp))
    print("TN :" + str(tn))
    print("FP :" + str(fp))
    print("FN :" + str(fn))