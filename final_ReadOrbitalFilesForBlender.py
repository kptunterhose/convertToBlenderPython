import os


TEST = True
BLENDER = False
if BLENDER:
    import bpy


COLOR = {
    'pos': (1, 0, 0, .25),
    'neg': (0, 0, 1, .25)
}

dictPosOrbitals = {}
dictNegOrbitals = {}
listOfNames = []

class Orbital(object):

    def __init__(self, filePath):
        self.role = []
        self.color = []
        self.vertices = []
        self.faces = []
        self.posOrbsVertices = []
        self.posOrbsFaces = []
        self.negOrbsVertices = []
        self.negOrbsFaces = []

        with open(filePath, 'r') as f:
            while True:
                line = f.readline()
                try:
                    firstWord = line.split()[0]
                except IndexError:
                    print('end of file reached')
                    break
                if firstWord == 'faces':
                    facesLine = f.readline().split()
                    if len(facesLine) == 3:
                        self.faces.append((int(facesLine[0]), int(facesLine[1]), int(facesLine[2])))
                    else:
                        print('something went wrong')
                elif firstWord == 'Verts':
                    threeLines = []
                    for i in range(3):
                        threeLines.append(f.readline())
                    # first line: color, role
                    self.color.append(int(threeLines[0].split()[1]))
                    self.role.append(int(threeLines[0].split()[3]))

                    # sedond line norm vector
                    # third line coords
                    coordsLine = threeLines[2].split()[1:]
                    self.vertices.append((float(coordsLine[0]), float(coordsLine[1]), float(coordsLine[2])))

                else:
                    break

    def makePosNegList(self):
        self.posList = []
        self.negList = []
        for i in range(len(self.role)):
            if self.role[i] == 0:
                self.posList.append(i)
            elif self.role[i] == 1:
                self.negList.append(i)
            else:
                print('something went wrong')

    def makePosOrbs(self):
        for i in self.posList:
            self.posOrbsVertices.append(self.vertices[i])
            for ii in range(len(self.faces)):
                if i in self.faces[ii]:
                    self.posOrbsFaces.append(self.faces[ii])

    def makeNegOrbs(self):
        try:
            smallesNeg = min(self.negList)
        except ValueError:
            print('no neg orbital exists')
        for i in self.negList:
            self.negOrbsVertices.append(self.vertices[i])
            for ii in range(len(self.faces)):
                if i in self.faces[ii]:
                    origFace = self.faces[ii]
                    newFace = (origFace[0] - smallesNeg, origFace[1] - smallesNeg, origFace[2] - smallesNeg)
                    self.negOrbsFaces.append(newFace)


if __name__ == '__main__':
    if TEST:
        basePath = 'orbitalTestData/'
    else:
        basePath = './'
    listOfNames = os.listdir(basePath)
    listOfOrbitalNames = []
    for name in listOfNames:
        if name.split('.')[-1] == 'orb':
            listOfOrbitalNames.append(name)

    for orbitalName in listOfOrbitalNames:
        print('reading orbital: ' + orbitalName)
        orb = Orbital(basePath + orbitalName)
        orb.makePosNegList()
        orb.makePosOrbs()
        orb.makeNegOrbs()
        dictPosOrbitals[orbitalName] = {}
        dictPosOrbitals[orbitalName]['verts'] = orb.posOrbsVertices
        dictPosOrbitals[orbitalName]['faces'] = orb.posOrbsFaces
        dictNegOrbitals[orbitalName] = {}
        dictNegOrbitals[orbitalName]['verts'] = orb.negOrbsVertices
        dictNegOrbitals[orbitalName]['faces'] = orb.negOrbsFaces

    print(len(dictPosOrbitals))
    print(len(dictNegOrbitals))

