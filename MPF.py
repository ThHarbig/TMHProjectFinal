from PDBHelixParser import PDBParser, PDBHelixParser
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from Helix import Helix
from matplotlib.colors import colorConverter
from ModuleMethods import  *
from mpl_toolkits.mplot3d import axes3d

__author__ = "Kevin Menden"
__date__ = '02.06.2016'
"""
This program takes a PDB file, reads it and predicts
which of its helices are transmembrane helices and which are not
Using this helices the program tries to find the most likely
orientation of the membrane and returns this membrane
as two planes: one for each ende of the membrane
"""


# Define the argparser



if __name__== "__main__":
    # parser = argparse.ArgumentParser(description="Membrane Plane Finder")
    # parser.add_argument('pdb')

    pdbParser = PDBHelixParser("test/1bm1.pdb")
    pdbParser.parse_pdb_file()
    print("Parsed PDB file succesfully")
    raw_helices = pdbParser.proteinHelixSequences

    # Convert raw helices into Helix objects
    helix_set = []
    for h in raw_helices:
        helix_set.append(Helix(h))
    print("Found " + str(len(helix_set)) + " helices")

    # Predict transmembrane helices
    tmh_set = []
    for h in helix_set:
        if is_transmembrane_helix(h) == 'tm':
            tmh_set.append(h)

    print("Predicted transmembrane helices: " + str(len(tmh_set)))

    # Calculate the normal vector of the transmembrane plane
    normal = calculate_normal_vector(tmh_set)

    # Choose the helix closest to the normal for width of the membrane
    # Create a plot
    lower_membrane = tmh_set[0].start_point
    upper_membrane = tmh_set[0].end_point

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for tmh in tmh_set:
        add_quiver(ax, tmh)

    ax.quiver(lower_membrane[0], lower_membrane[1], lower_membrane[2], normal[0], normal[1], normal[2],
              colors=colorConverter.to_rgb('g'), lw=8, length=50, pivot='tail')
    d1 = -lower_membrane.dot(normal)
    d2 = -upper_membrane.dot(normal)
    xx1, yy1 = np.meshgrid(range(100), range(100))
    xx2, yy2 = np.meshgrid(range(100), range(100))
    z1 = (-normal[0] * xx1 - normal[1] * yy1 - d1) * 1. / normal[2]
    z2 = (-normal[0] * xx2 - normal[1] * yy2 - d2) * 1. / normal[2]
    ax.plot_surface(xx1, yy1, z1, cmap='BuPu', lw=4)
    ax.plot_surface(xx2, yy2, z2, cmap='Blues', lw=4)

    ax.set_xlim(-50, 50)
    ax.set_ylim(-50, 50)
    ax.set_zlim(-50, 50)

    plt.show()

