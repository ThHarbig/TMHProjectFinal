import numpy as np



class MembranePlacer(object):
    def __init__(self, helices, structure, normal, file):
        self.normal = normal
        self.helices = helices
        self.helicalPos = []
        self.structure = structure
        self.file = file
        self.stepsize=np.array([0,0,0])
        self.score = 0
        self.tmhPos=[]
        for helix in self.helices:
            if self.closeToNormal(helix.vector):
                self.tmhPos.extend(helix.positions)

    def computeMiddlePlane(self):
        averagemiddle = np.array([0., 0., 0.])
        counter=0
        for helix in self.helices:
            start = helix.start_point
            end = helix.end_point
            if self.closeToNormal(helix.vector):
                if helix.length < 50:
                    averagemiddle[0] += ((start[0] + end[0]) / 2)
                    averagemiddle[1] += ((start[1] + end[1]) / 2)
                    averagemiddle[2] += ((start[2] + end[2]) / 2)
                    counter+=1
        for i in range(0, len(averagemiddle)):
            averagemiddle[i] /= counter
        return averagemiddle

    def findMembrane(self):
        up=self.shiftUp()
        down=self.shiftDown()
        maximum=np.max([up[2],down[2]])
        if up[2]==maximum:
            return up[0],up[1]
        else:
            return down[0],down[1]

    def shiftMembrane(self,plane,distance):
        return plane+distance

    def shiftUp(self):
        middle=self.computeMiddlePlane()
        A=self.shiftMembrane(middle,15*self.stepsize)
        B=self.shiftMembrane(middle,-15*self.stepsize)
        x1, x2, x3, y1, y2, y3, score = 0,0,0,0,0,0,0
        for i in range(15):
            A=self.shiftMembrane(A,-self.stepsize)
            B =self.shiftMembrane(B,-self.stepsize)
            broadened=self.broaden(A,B)
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
        x1, x2, x3, y1, y2, y3, score = 0,0,0,0,0,0,0
        for i in range(15):
            middle = self.computeMiddlePlane()
            A = self.shiftMembrane(middle, 15 * self.stepsize)
            B = self.shiftMembrane(middle, -15 * self.stepsize)
            broadened=self.broaden(A, B)
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
        x1, x2, x3, y1, y2, y3, score = 0,0,0,0,0,0,0
        for i in range(10):
            A=self.shiftMembrane(A,self.stepsize)
            B =self.shiftMembrane(B,-self.stepsize)
            newScore=self.scoring(A,B)
            if newScore>score:
                score=newScore
                x1=A[0]
                x2=A[1]
                x3=A[2]
                y1=B[0]
                y2=B[1]
                y3=B[2]
        return np.array([x1,x2,x3]),np.array([y1,y2,y3]),score


    def helicalPositions(self):
        pdb = open(self.file)
        lines = pdb.readlines()
        for line in lines:
            if line.startswith("HELIX"):
                if line[39:40] == "1":
                    self.helicalPos.extend(range(int(line[22:25]), int(line[34:37])))

    def computeStepSize(self):
        self.stepsize= self.normalize(self.normal)

    def normalize(self,vector):
        length=np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
        return vector/length

    def closeToNormal(self,vec1):
        angle=(np.arccos(self.normalize(vec1).dot(self.normalize(self.normal))))
        if not(angle>np.pi/3. and angle<2*(np.pi/3.) or (angle>4*(np.pi/3.) and angle<5*(np.pi/3.))):
            return True
        else:
            return False

    def scoring(self, A, B):
        hydrophobic_residues = ["PHE", "GLY", "ILE", "LEU", "MET", "VAL", "TRP", "THR"]
        hydrophilic_residues = ["ALA", "CYS", "ASP", "GLU", "HIS", "LYS", "ASN", "PRO", "GLN", "ARN", "SER", "THR"]
        hphobcount = 0
        hphilcount = 0
        helixcount = 0
        loopcount = 0
        tmhcount=0
        notTmhCount=0
        for residue in self.structure.get_residues():
            try:
                atom = residue['CA']
                coordinates = atom.get_coord()
                up = coordinates - B
                down = coordinates - A
                first = up.dot(B)
                second = down.dot(A)
                if (first < 0 and second > 0) or (first > 0 and second < 0):
                    if residue.get_resname() in hydrophilic_residues:
                        if residue.id[1] in self.helicalPos:
                            hphilcount += 1
                    if residue.get_resname in hydrophilic_residues:
                        if residue.id[1] in self.helicalPos:
                            hphobcount += 1
                    if residue.id[1] in self.helicalPos:
                        helixcount += 1
                    else:
                        loopcount += 1
                    if residue.id[1] in self.tmhPos:
                        tmhcount+=1
                    else:
                        notTmhCount+=1

            except KeyError:
                continue

        return (hphobcount  / (hphilcount + 1))+4*(tmhcount/(notTmhCount+1)+(helixcount)/(loopcount+1))
        #return ((helixcount)/(loopcount+1))
        #return (tmhcount/(notTmhCount+1))
        #return ((hphobcount ) / (hphilcount + 1))

    def positiveInsideRule(self,A,B):
        positiveAS=["ARG","LYS","HIS"]
        side1pos=0
        side1nonpos=0
        side2pos=0
        side2nonpos=0
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
                        side1pos+=1
                    else:
                        side1nonpos+=1
                if first > 0 and second > 0:
                    if residue.get_resname() in positiveAS:
                        side2pos += 1
                    else:
                        side2nonpos += 1
            except KeyError:
                continue
        if side1pos/side1nonpos+1>side2pos/side2nonpos+1:
            return 0
        else:
            return 1

    def placeMembrane(self):
        self.helicalPositions()
        self.computeStepSize()
        result=self.findMembrane()
        positive=self.positiveInsideRule(result[0],result[1])
        return result,positive
