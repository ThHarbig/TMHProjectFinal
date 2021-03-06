from PDBHelixParser import PDBParser, PDBHelixParser
import matplotlib.pyplot as plt
from Helix import Helix
from matplotlib.colors import colorConverter
from ModuleMethods import  *
from Membrane import MembranePlacer
import warnings
import argparse
import sys
from Bio import BiopythonWarning
from mpl_toolkits.mplot3d import axes3d

__date__ = '02.06.2016'
"""
This program takes a PDB file, reads it and predicts
which of its helices are transmembrane helices and which are not
Using this helices the program tries to find the most likely
orientation of the membrane and returns this membrane
as two planes: one for each end of the membrane
"""

# Define the argparser
if __name__== "__main__":
    parser = argparse.ArgumentParser(description="Prediction of Transmembrane Helices")
    parser = argparse.ArgumentParser(description="Membrane Plane Finder")
    parser.add_argument("pdb_file", help="choose pdb file to predict", nargs='?')
    parser.add_argument("plot", help="toggles the plot", nargs='?', default=False, type=bool)
    args = parser.parse_args()
    if args.pdb_file is None:
        print("Missing pdb file input!")
        sys.exit()
    file = args.pdb_file
    #ignore warnings regarding discontinuous chainsr
    warnings.simplefilter('ignore', BiopythonWarning)
    # parse the given pdb file
    pdbParser = PDBHelixParser(file)
    pdbParser.parse_pdb_file()
    structure = pdbParser.structure         # The whole structure of the PDB file
    print("Parsed PDB file succesfully.")
    raw_helices = pdbParser.proteinHelixSequences
    # Convert raw helices into Helix objects
    helix_set = []
    for h in raw_helices:
        if len(h) > 0:
            helix_set.append(Helix(h))
    print("Found " + str(len(helix_set)) + " helices.")

    # Predict transmembrane helices
    tmh_set = []
    for h in helix_set:
        if is_transmembrane_helix(h) == 'tm':
            tmh_set.append(h)
    #print(len(tmh_set) / len(helix_set))
    if len(tmh_set) < 3 or len(tmh_set)/len(helix_set) < 0.1:
        print("The given protein is no transmembrane protein!")
        sys.exit(0)
    else:
        print("The given protein is a transmembrane protein!")

    print("Predicted transmembrane helices: " + str(len(tmh_set)))

    # Calculate the normal vector of the transmembrane plane
    normal = calculate_normal_vector(tmh_set)

    # Choose the helix closest to the normal for width of the membrane
    # Create a plot
    placer=MembranePlacer(tmh_set,structure,normal,file)
    membranes=placer.placeMembrane()
    lower_membrane = membranes[0]
    upper_membrane = membranes[1]
    thickness = ((normal.dot(upper_membrane) - normal.dot(lower_membrane))) / 2
    if thickness < 0:
        thickness = thickness*-1
    print("Normal vector of the transmembrane plane: " + str(normal))
    print("Thicknes of the membrane: " + str(thickness) + " A")
    if not placer.proportionOfHelices(lower_membrane,upper_membrane):
        print("WARNING: Membrane could possibly be not placed correctly. More than 40% of the transmembrane helices are outside of the membrane planes")
    print("Position of the first membrane: " + str(lower_membrane))
    print("Position of the second membrane: " + str(upper_membrane))


    if(args.plot):
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
        ax.plot_surface(xx1, yy1, z1, color="red", lw=4, )
        ax.plot_surface(xx2, yy2, z2, color="red", lw=4)
        ax.set_xlim(-100, 100)
        ax.set_ylim(-100, 100)
        ax.set_zlim(-100, 100)
        plt.show()
