from xml.dom.minidom import parse
from Bio.PDB import *


def parseXML(file):
    sequences = []
    data = parse(file)
    all = data.getElementsByTagName("pdbtm")
    count = 0
    for pdbtm in all:
        if count == 230:
            break
        chains = pdbtm.getElementsByTagName("CHAIN")
        for chain in chains:
            if chain.getAttribute("TYPE") == "alpha":
                count += 1
                sequences.extend(searchHelices(pdbtm))
                break
    return sequences


def searchHelices(pdbtm):
    results = []
    try:
        pdbHelices = getHelices(pdbtm.getAttribute("ID"))
        chains = pdbtm.getElementsByTagName("CHAIN")
        for chain in chains:
            if chain.getAttribute("TYPE") == "alpha":
                id = chain.getAttribute("CHAINID")
                sequence = chain.getElementsByTagName("SEQ")[0].firstChild.nodeValue.replace(" ", "").replace("\n", "")
                regions = chain.getElementsByTagName("REGION")
                for region in regions:
                    if region.getAttribute("type") == "H":
                        start = int(region.getAttribute("pdb_beg"))
                        if id in pdbHelices:
                            for entry in pdbHelices[id]:
                                if start >= entry["start"]:
                                    if start <= entry["end"]:
                                        subseq = sequence[entry["start"]:entry["end"]]
                                        if "U" not in subseq:
                                            results.append(subseq)
                                        break
    except IOError:
        print("file not found")
    return results


def getHelices(file):
    dict = {}
    pdb = open("tmh_set/pdb" + file + ".ent", "r")
    lines = pdb.readlines()
    for line in lines:
        if line.startswith("HELIX"):
            if line[39:40] == "1":
                chain = line[19]
                if not chain in dict:
                    dict[chain] = []
                dict[chain].append({"start": int(line[22:25]), "end": int(line[34:37])})
    return dict


if __name__ == "__main__":
    sequences = parseXML("pdbtmall.xml")
    print(len(sequences))