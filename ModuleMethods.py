import math
import numpy as np
from scipy.spatial import distance
__author__ = "Kevin Menden"
__date__ = '01.06.2016'

"""
Methods
"""

from sklearn.decomposition import PCA
from sklearn import svm
from sklearn.externals import joblib

hydrophobic_residues = ["F", "G", "I", "L", "M", "V", "W", "Y"]

def calculate_normal_vector(helices):
    """
    Calculate the mean vector of all helix vectors, i.e. the
    membrane normal vector
    :param helices: all membrane helices
    :return: the normal vector
    """
    helix_set = []
    for h in helices:
        helix_set.append(h.vector)
        helix_set.append(h.neg_vector)
    pca = PCA(n_components=1)
    pca.fit_transform(helix_set)
    return pca.components_[0]

def add_quiver(ax, helix):
    """
    Add a quiver to the plot
    :param ax:
    :param helix:
    :return:
    """
    vec = helix.vector
    start = helix.start_point
    end = helix.end_point
    dist = end - start
    length = math.sqrt(dist[0]**2 + dist[1]**2 + dist[2]**2)
    ax.quiver(start[0], start[1], start[2], vec[0], vec[1], vec[2], length=length, arrow_length_ratio=0.1, pivot='tail',
              cmap='Accent', lw=3)

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

def is_transmembrane_helix(helix):
    """
    Predicts if a helix is transmembrane or not
    :param helix: the helix of interest (sequence)
    :return: True or False
    """
    clf = joblib.load("tmh_predictor_weighted.pkl")
    feats = calculate_features(helix.sequence)
    return clf.predict([feats])

def normalize(vector):
    """
    normalizes a vector
    :param vector:
    :return:
    """
    length = np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
    return vector / length


