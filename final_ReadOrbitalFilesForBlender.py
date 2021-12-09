"""
Copy this script in the folder with the .orb fils.
If you don't want all orbitals add option -o and after its the list of wanted orbitals
e.g. python final_ReadOrbitalFilesForBlender.py -o 9 12 15

"""
#TODO: check ob orbitals n surfaces / vorzeichen -> entsprechende objekte trennen -> für modifikation


import sys
import os
import json
import queue as q
TEST = True

listOfNames = []
dOrbs = {}
content1 = '''
import bpy

CollectionName = 'Orbital'

dOrbs = '''

content2 = """
material_pos = bpy.data.materials.new(name='material-pos')
material_pos.diffuse_color = (1,0,0,.25)
material_neg = bpy.data.materials.new(name='material-neg')
material_neg.diffuse_color = (0,0,1,.25)

material = {
    'pos': material_pos,
    'neg': material_neg
}


try: 
    bpy.data.collections['Orbital']
except KeyError:
    orbitalCollection = bpy.data.collections.new(CollectionName)
    bpy.context.scene.collection.children.link(orbitalCollection)


for orbital in dOrbs:
    try:
        bpy.data.collections[str(orbital)]
    except KeyError:
        singleOrbital = bpy.data.collections.new(str(orbital))
        bpy.context.scene.collection.children[CollectionName].children.link(singleOrbital)
    for sign in dOrbs[orbital]:
        dOrbs[orbital][sign]['mesh'] = bpy.data.meshes.new(str(orbital) + sign + '_mesh')
        dOrbs[orbital][sign]['mesh'].from_pydata(dOrbs[orbital][sign]['verts'], [], dOrbs[orbital][sign]['faces'])
        dOrbs[orbital][sign]['obj'] = bpy.data.objects.new(sign , dOrbs[orbital][sign]['mesh'])
        dOrbs[orbital][sign]['obj'].data.materials.append(material[sign.split('-')[0]])
        dOrbs[orbital][sign]['obj'].scale = [.5, .5, .5]
        singleOrbital.objects.link(dOrbs[orbital][sign]['obj'])
"""


def readOrbitalFromFiles(listOfOrbitalNames, onlySomeOrbitals = []):
    for orbitalName in listOfOrbitalNames:
        add = False
        if len(onlySomeOrbitals) == 0:
            add = True
        else:
            for onlySome in onlySomeOrbitals:
                if 'o'+ onlySome in orbitalName:
                    add = True
        if add:
            dOrbs[orbitalName] = {}

            print('reading orbital: ' + orbitalName)

            orb = Orbital(basePath + orbitalName)
            orb.makePosNegList()
            orb.makePosOrbs()
            orb.makeNegOrbs()

            orb.doFacesStuff()
            dictOfFacesPos, dictOfVertsPos = orb.getPostiveFacesAndVertsDict()
            dictOfFacesNeg, dictOfVertsNeg = orb.getNegativeFacesAndVertsDict()

            if not orb.testStuff():
                exit('something messd up with faces and verts')

            for nSurface in dictOfFacesPos:
                dOrbs[orbitalName]['pos-' + str(nSurface)] = {
                    'verts': dictOfVertsPos[nSurface],
                    'faces': dictOfFacesPos[nSurface]
                }
            for nSurface in dictOfFacesNeg:
                dOrbs[orbitalName]['neg-' + str(nSurface)] = {
                    'verts': dictOfVertsNeg[nSurface],
                    'faces': dictOfFacesNeg[nSurface]
                }


def makePythonOutfile(outFileName):
    outFile = open(outFileName, 'w')
    outFile.write(content1)
    # write orbital stuff here
    outFile.write(json.dumps(dOrbs))
    outFile.write(content2)
    outFile.close()
    print('write output to ' + outFileName)
    return 0


def getBesucht(adj, start, suche):
    # adj ist die Adjazenzliste {knoten: [kanten]}
    # start ist der Index des Knoten, in dem die Suche beginnt
    # suche ist der gesuchte Knoten
    queue = q.Queue()
    queue.put(start)
    besucht = set()
    while queue.qsize() > 0:
        aktiverKnoten = queue.get()
        print(aktiverKnoten)
        #if aktiverKnoten in besucht:
        #    pass
        #else:
        besucht.add(aktiverKnoten)
        for andererKnoten in adj[aktiverKnoten]:
            if andererKnoten in besucht:
                continue
            if andererKnoten == suche:
                # Knoten gefunden
                return besucht
            queue.put(andererKnoten)
    print(besucht)
    return besucht


def getBesuchtNeu(adjazenzlist, start):
    queue = q.Queue()
    queue.put(start)
    besucht = set()
    while queue.qsize() > 0:
        aktiverKnoten = queue.get()
        besucht.add(aktiverKnoten)
        for andererKnoten in adjazenzlist[aktiverKnoten]:
            if andererKnoten not in queue.queue and andererKnoten not in besucht:
                queue.put(andererKnoten)
    return besucht


def getDictOfKnots(listOfFaces, listOfVerts):
    knotenDict = dict()
    for i in range(len(listOfVerts)):
        knotenDict[i] = set()
        for triangle in listOfFaces:
            if i in triangle:
                for knot in triangle:
                    if knot != i:
                        knotenDict[i].add(knot)
    return knotenDict


def getDictOfSurfaces(knoten):
    listOfKnots = list(knoten.keys())
    surfaces = {}
    nSurface = 0
    while listOfKnots:
        #s = getBesucht(knoten, listOfKnots[0], -1)
        s = getBesuchtNeu(knoten, listOfKnots[0])
        for k in s:
            listOfKnots.remove(k)
        surfaces[nSurface] = s
        nSurface += 1
    return surfaces


class Orbital(object):

    def __init__(self, filePath):
        self.role = list()
        self.color = list()
        self.vertices = list()
        self.faces = list()
        self.posOrbsVertices = list()
        self.posOrbsFaces = list()
        self.negOrbsVertices = list()
        self.negOrbsFaces = list()
        self.dKnotsPos = dict()
        self.dKnotsNeg = dict()
        self.dictOfFacesPos = dict()
        self.dictOfFacesNeg = dict()
        self.dictOfVertsPos = dict()
        self.dictOfVertsNeg = dict()
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

    def makeDictOfFacesAndVerts(self, faces, verts, dictOfSurfaces):
        dictOfFaces = dict()
        dictOfVerts = dict()
        mapOfIndices = dict()  # {oldIndex: newIndex}
        for surface in dictOfSurfaces:
            dictOfFaces[surface] = []
            for face in faces:
                for point in face:
                    if point in dictOfSurfaces[surface]:
                        dictOfFaces[surface].append(list(face))
                        break
            dictOfVerts[surface] = []
            for knot in dictOfSurfaces[surface]:
                dictOfVerts[surface].append(verts[knot])

            # mapping of indices
            indexNew = 0
            actualSurface = list(dictOfSurfaces[surface])  # entries in a set are not sorted
            actualSurface.sort()
            for indexOld in actualSurface:
                mapOfIndices[indexOld] = indexNew
                indexNew += 1

            for i in range(len(dictOfFaces[surface])):
                for ii in range(len(dictOfFaces[surface][i])):
                    dictOfFaces[surface][i][ii] = mapOfIndices[dictOfFaces[surface][i][ii]]

        return dictOfFaces, dictOfVerts

    def makeDictOfPosNegFacesAndVerts(self):
        self.dictOfFacesPos, self.dictOfVertsPos = self.makeDictOfFacesAndVerts(self.posOrbsFaces,
                                                                           self.posOrbsVertices,
                                                                           self.dSurfacePos)
        self.dictOfFacesNeg, self.dictOfVertsNeg = self.makeDictOfFacesAndVerts(self.negOrbsFaces,
                                                                           self.negOrbsVertices,
                                                                           self.dSurfaceNeg)
        pass

    def makeDictOfSurfaces(self, knoten=dict()):
        # takes long runtime -> optimize !!!
        listOfKnots = list(knoten.keys())
        surfaces = {}
        nSurface = 0
        while listOfKnots:
            #s = getBesucht(knoten, listOfKnots[0], -1)
            s = getBesuchtNeu(knoten, listOfKnots[0])
            for k in s:
                listOfKnots.remove(k)
            surfaces[nSurface] = s
            nSurface += 1
        return surfaces

    def makeDictOfPosNegSurfaces(self):
        self.dSurfacePos = self.makeDictOfSurfaces(self.dKnotsPos)
        self.dSurfaceNeg = self.makeDictOfSurfaces(self.dKnotsNeg)
        pass

    def makeDictOfKnots(self, listOfFaces, listOfVerts):
        knotenDict = dict()
        for i in range(len(listOfVerts)):
            knotenDict[i] = set()
            for triangle in listOfFaces:
                if i in triangle:
                    for knot in triangle:
                        if knot != i:
                            knotenDict[i].add(knot)
        return knotenDict

    def makeDictOfPosNegKnots(self):
        self.dKnotsPos = self.makeDictOfKnots(self.posOrbsFaces, self.posOrbsVertices)
        self.dKnotsNeg = self.makeDictOfKnots(self.negOrbsFaces, self.negOrbsVertices)
        pass

    def doFacesStuff(self):
        self.makeDictOfPosNegKnots()
        self.makeDictOfPosNegSurfaces()
        self.makeDictOfPosNegFacesAndVerts()
        pass

    def getPostiveFacesAndVertsDict(self):
        return self.dictOfFacesPos, self.dictOfVertsPos

    def getNegativeFacesAndVertsDict(self):
        return self.dictOfFacesNeg, self.dictOfVertsNeg

    def testStuff(self):
        #check Anzahl Vertices mit max Wert von Faces
        #get Anzahl vertices
        result = True
        nVertsPos = dict()
        for s in self.dictOfVertsPos:
            nVertsPos[s] = len(self.dictOfVertsPos[s])
        nVertsNeg = dict()
        for s in self.dictOfVertsNeg:
            nVertsNeg[s] = len(self.dictOfVertsNeg[s])
        #get max Wert von Faces
        maxVertPos = dict()
        for s in self.dictOfFacesPos:
            for p in self.dictOfFacesPos[s]:
                maxVertPos[s] = max(p) + 1
        maxVertNeg = dict()
        for s in self.dictOfFacesNeg:
            for p in self.dictOfFacesNeg[s]:
                maxVertNeg[s] = max(p) + 1

        for s in nVertsPos:
            result &= (nVertsPos[s] == maxVertPos[s])
        for s in nVertsNeg:
            result &= (nVertsNeg[s] == maxVertNeg[s])
        return result


if __name__ == '__main__':
    onlyList = []
    if len(sys.argv) > 2:
        if '-o' in sys.argv:
            onlyList = sys.argv[sys.argv.index('-o')+1:]
    if onlyList:
        print('making following orbitals', onlyList)
    else:
        print('making all orbitals i could find')
    if TEST:
        basePath = 'orbitalTestData/'
    else:
        basePath = './'
    listOfNames = os.listdir(basePath)
    listOfOrbitalNames = []
    for name in listOfNames:
        if name.split('.')[-1] == 'orb':
            listOfOrbitalNames.append(name)

    readOrbitalFromFiles(listOfOrbitalNames, onlyList)
    makePythonOutfile(basePath + 'drawOrbitals_blender.py')


