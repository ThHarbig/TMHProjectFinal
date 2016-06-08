import numpy as np


class MembranePlacer(object):
    def __init__(self, helices, structure, normal, file):
        self.normal = normal
        self.helices = helices
        self.helicalPos = []
        self.structure = structure
        self.file = file
        self.score = 0
        self.tmhPos = []
        for helix in self.helices:
            if self.closeToNormal(helix.vector):
                if helix.length < 50 and helix.length > 15:
                    self.tmhPos.extend(helix.positions)

    def computeMiddlePlane(self):
        """
        Computes a plane using the average midpoint of all helices
        (helicies with a "wrong" angle and helices which are too long are excluded)
        :return: middle plane
        """
        averagemiddle = np.array([0., 0., 0.])
        counter = 0
        for helix in self.helices:
            start = helix.start_point
            end = helix.end_point
            #the angle of the helix to the normal is too close to 90°
            if self.closeToNormal(helix.vector):
                #helix too long or too short
                if helix.length < 50 and helix.length > 15:
                    averagemiddle[0] += ((start[0] + end[0]) / 2)
                    averagemiddle[1] += ((start[1] + end[1]) / 2)
                    averagemiddle[2] += ((start[2] + end[2]) / 2)
                    counter += 1
        for i in range(0, len(averagemiddle)):
            averagemiddle[i] /= counter
        return averagemiddle

    def findMembrane(self):
        """
        Finds the best positions for the membrane by shifting the two planes up and down
        :return: best membrane planes
        """
        up = self.shiftUp()
        down = self.shiftDown()
        maximum = np.max([up[2], down[2]])
        if up[2] == maximum:
            return up[0], up[1]
        else:
            return down[0], down[1]

    def shiftUp(self):
        """
        Computes the first guess for the two planes, shifts them up and broadens them in each step
        :return: best membrane planes (by shifting up)
        """
        middle = self.computeMiddlePlane()
        A = middle + 15 * self.normal
        B = middle - 15 * self.normal
        x1, x2, x3, y1, y2, y3, score = 0, 0, 0, 0, 0, 0, 0
        for i in range(0, 15):
            A = A + self.normal
            B = B + self.normal
            broadened = self.broaden(A, B)
            if broadened[2] > score:
                score = broadened[2]
                x1 = broadened[0][0]
                x2 = broadened[0][1]
                x3 = broadened[0][2]
                y1 = broadened[1][0]
                y2 = broadened[1][1]
                y3 = broadened[1][2]
        return np.array([x1, x2, x3]), np.array([y1, y2, y3]), score

    def shiftDown(self):
        """
        Computes the first guess for the two planes, shifts them down and broadens them in each step
        :return: best membrane planes (by shifting down)
        """
        middle = self.computeMiddlePlane()
        A = middle + 15 * self.normal
        B = middle - 15 * self.normal
        x1, x2, x3, y1, y2, y3, score = 0, 0, 0, 0, 0, 0, 0
        for i in range(0, 15):
            A = A - self.normal
            B = B - self.normal
            broadened = self.broaden(A, B)
            if broadened[2] > score:
                score = broadened[2]
                x1 = broadened[0][0]
                x2 = broadened[0][1]
                x3 = broadened[0][2]
                y1 = broadened[1][0]
                y2 = broadened[1][1]
                y3 = broadened[1][2]
        return np.array([x1, x2, x3]), np.array([y1, y2, y3]), score

    def broaden(self, A, B):
        """
        increases the distance of the two membranes
        :param A:membrane plane
        :param B:membrane plane
        :return: best membrane (by broadening)
        """
        x1, x2, x3, y1, y2, y3, score = 0, 0, 0, 0, 0, 0, 0
        for i in range(10):
            A = A + self.normal
            B = B - self.normal
            newScore = self.scoring(A, B)
            if newScore > score:
                score = newScore
                x1 = A[0]
                x2 = A[1]
                x3 = A[2]
                y1 = B[0]
                y2 = B[1]
                y3 = B[2]
        return np.array([x1, x2, x3]), np.array([y1, y2, y3]), score

    def helicalPositions(self):
        """
        Gets all the positions of the protein which are in a helix
        :return:
        """
        pdb = open(self.file)
        lines = pdb.readlines()
        for line in lines:
            if line.startswith("HELIX"):
                if line[39:40] == "1":
                    self.helicalPos.extend(range(int(line[22:25]), int(line[34:37])))

    def normalize(self, vector):
        """
        normalizes a vector
        :param vector:
        :return:
        """
        length = np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
        return vector / length

    def closeToNormal(self, vec1):
        """
        Checks if a vector is close to the normal
        :param vec1:
        :return:
        """
        angle = (360 / (2 * np.pi)) * (np.arccos(self.normalize(vec1).dot(self.normal)))
        if angle < 60 or angle > 120:
            return True
        else:
            return False

    def scoring(self, A, B):
        """
        scores the positions of the two membranes
        :param A:Membrane plane
        :param B:Membrane plane
        :return:score
        """
        hydrophobic_residues = ["PHE", "GLY", "ILE", "LEU", "MET", "VAL", "TRP", "THR"]
        hydrophilic_residues = ["ALA", "CYS", "ASP", "GLU", "HIS", "LYS", "ASN", "PRO", "GLN", "ARN", "SER", "THR"]
        hphobcount, hphilcount, helixcount, loopcount, tmhcount, notTmhCount = 0, 0, 0, 0, 0, 0
        #loop through all residues
        for residue in self.structure.get_residues():
            try:
                atom = residue['CA']
                coordinates = atom.get_coord()
                distA = (self.normal.dot(coordinates) - self.normal.dot(A)) / (self.normal.dot(self.normal))
                distB = (self.normal.dot(coordinates) - self.normal.dot(B)) / (self.normal.dot(self.normal))
                #check if residue is between the planes
                if (distB < 0 < distA) or (distB > 0 > distA):
                    if residue.get_resname() in hydrophobic_residues:
                        hphobcount += 1
                    if residue.get_resname() in hydrophilic_residues:
                        hphilcount += 1
                    if residue.id[1] in self.helicalPos:
                        helixcount += 1
                    else:
                        loopcount += 1
                    if residue.id[1] in self.tmhPos:
                        tmhcount += 1
                    else:
                        notTmhCount += 1
            except KeyError:
                continue
        #scoring function:
        #return (hphobcount/((hphilcount+1)))
        return ((hphobcount) / (hphilcount + 1)) + ((tmhcount) / (notTmhCount + 1)+(helixcount/(loopcount+1)))

    def positiveInsideRule(self, A, B):
        positiveAS = ["ARG", "LYS", "HIS"]
        side1pos = 0
        side1nonpos = 0
        side2pos = 0
        side2nonpos = 0
        for residue in self.structure.get_residues():
            try:
                atom = residue['CA']
                coordinates = atom.get_coord()
                up = coordinates - A
                down = coordinates - B
                first = up.dot(A)
                second = down.dot(B)
                if first < 0 and second < 0:
                    if residue.get_resname() in positiveAS:
                        side1pos += 1
                    else:
                        side1nonpos += 1
                if first > 0 and second > 0:
                    if residue.get_resname() in positiveAS:
                        side2pos += 1
                    else:
                        side2nonpos += 1
            except KeyError:
                continue
        if side1pos / (side1nonpos + 1) > side2pos / (side2nonpos + 1):
            return 0
        else:
            return 1

    def placeMembrane(self):
        self.helicalPositions()
        result = self.findMembrane()
        return result
