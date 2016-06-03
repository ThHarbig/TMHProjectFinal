from PDBextractor import SecondaryStructureCounter

__author__ = "Kevin Menden"
__date__ = '24.05.2016'

from sklearn import svm, cross_validation
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB.PDBParser import PDBParser
from XMLParser import parseXML
from sklearn.externals import joblib
from matplotlib import cm as cm



hydrophobic_residues = ["F", "G", "I", "L", "M", "V", "W", "Y"]
hydrophilic_residues = ["A", "C", "D", "E", "H", "K", "N", "P", "Q", "R", "S", "T"]

def calculate_hydrophobicity(seq):
    """
    Calculate the hydrophobicity factor for an amino acid sequence
    :param seq: the input amino acid sequence
    :return: the hydrophobicity factor
    """
    phobic = 0
    seq = seq.upper()
    for elem in seq:
        if elem in hydrophobic_residues:
            phobic += 1

    # if len(seq) > 0:
    #
    #     return phobic/len(seq)
    # else: return 0
    return phobic




def calculate_features(sequence):
    """
    Calculates all features of a given amino acid sequence
    :param sequence: amino acid sequence
    :return: array of features
    """
    hp_factor = calculate_hydrophobicity(sequence)
    length_factor = len(sequence)

    return [hp_factor, length_factor]

def make_training_set(tm_set, glob_set):
    """
    Make a training set in correct format for the SVM
    :param tm_set:
    :param glob_set:
    :return: the training set, the set of class labels
    """
    training = []
    lables = []

    # Append transmembrane set
    for elem in tm_set:
        training.append(calculate_features(elem))
        lables.append("tm")
    # Appen globular set
    for elem in glob_set:
        training.append(calculate_features(elem))
        lables.append("glob")

    return training, lables



tm_set = []

glob_set = []

# Parse PDB tmh_set and extract helices
structureCounter = SecondaryStructureCounter("globular_set")
structureCounter.parse_all_pdb_files()

new_helices = structureCounter.proteinHelixSequences
for prot in new_helices:
    for res in prot:
        glob_set.append(structureCounter.convertResiduesToLetters(res))
        # print(structureCounter.convertResiduesToLetters(res))

# Parse XML file containing helices
alpha_helices = parseXML("pdbtmall.xml")
print(len(alpha_helices))
for i in range(2500):
    tm_set.append(alpha_helices[i])


training, class_labels = make_training_set(tm_set, glob_set)

# Training
clf = svm.SVC(kernel='rbf')
clf.fit(X=training, y=class_labels)

# Cross-validation
print("Accuracy: " + str(np.mean(cross_validation.cross_val_score(clf, X=training, y=class_labels))))


# Make plot
x1 = []
y1 = []
x2 = []
y2 = []
for elem in tm_set:
    tmp = calculate_features(elem)
    x1.append(tmp[0])
    y1.append(tmp[1])

for elem in glob_set:
    tmp = calculate_features(elem)
    x2.append(tmp[0])
    y2.append(tmp[1])

fig = plt.figure()

x1 = np.array(x1)
y1 = np.array(y1)
ax1 = fig.add_subplot(211)
ax1.hexbin(x1, y1, gridsize=100, bins=None, cmap=cm.jet)
# plt.hexbin(x1, y1, gridsize=100, bins=None, cmap=cm.jet)
ax1.axis([x1.min(), x1.max(), y1.min(), y1.max()])
# ax1.colorbar()

ax2 = fig.add_subplot(212)
x2 = np.array(x2)
y2 = np.array(y2)
ax2.hexbin(x2, y2, gridsize=50, bins=None, cmap=cm.jet)
ax2.axis([x1.min(), x1.max(), y1.min(), y1.max()])

# ax2.colorbar()
plt.show()


# plt.xlabel("Hydrophobicity factor")
# plt.ylabel("Length")
# plt.show()

# Save the model
# joblib.dump(clf, "tmh_predictor.pkl")


