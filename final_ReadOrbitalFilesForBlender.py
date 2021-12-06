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
            dKnotsPos = getDictOfKnots(orb.posOrbsFaces, orb.posOrbsVertices)
            dKnotsNeg = getDictOfKnots(orb.negOrbsFaces, orb.negOrbsVertices)

            dictOfSurfacesPos = getDictOfSurfaces(dKnotsPos)
            dictOfSurfacesNeg = getDictOfSurfaces(dKnotsNeg)

            dictOfFacesPos = {}
            dictOfFacesNeg = {}
            dictOfVertsPos = {}
            dictOfVertsNeg = {}

            for surface in dictOfSurfacesPos:
                dictOfFacesPos[surface] = []
                for face in orb.posOrbsFaces:
                    for point in face:
                        if point in dictOfSurfacesPos[surface]:
                            dictOfFacesPos[surface].append(list(face))
                            break
                dictOfVertsPos[surface] = []
                for knot in dictOfSurfacesPos[surface]:
                    dictOfVertsPos[surface].append(orb.posOrbsVertices[knot])
                for i in range(len(dictOfFacesPos[surface])):
                    for ii in range(len(dictOfFacesPos[surface][i])):
                        dictOfFacesPos[surface][i][ii] -= min(dictOfSurfacesPos[surface])

            for surface in dictOfSurfacesNeg:
                dictOfFacesNeg[surface] = []
                for face in orb.negOrbsFaces:
                    for point in face:
                        if point in dictOfSurfacesNeg[surface]:
                            dictOfFacesNeg[surface].append(list(face))
                            break
                dictOfVertsNeg[surface] = []
                for knot in dictOfSurfacesNeg[surface]:
                    dictOfVertsNeg[surface].append(orb.negOrbsVertices[knot])
                for i in range(len(dictOfFacesNeg[surface])):
                    for ii in range(len(dictOfFacesNeg[surface][i])):
                        dictOfFacesNeg[surface][i][ii] -= min(dictOfSurfacesNeg[surface])

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
        besucht.add(aktiverKnoten)
        for andererKnoten in adj[aktiverKnoten]:
            if andererKnoten in besucht:
                continue
            if andererKnoten == suche:
                # Knoten gefunden
                return besucht
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
        s = getBesucht(knoten, listOfKnots[0], -1)
        for k in s:
            listOfKnots.remove(k)
        surfaces[nSurface] = s
        nSurface += 1
    return surfaces


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
    onlyList = []
    if len(sys.argv) > 2:
        if '-o' in sys.argv:
            onlyList = sys.argv[sys.argv.index('-o')+1:]
    print(onlyList)
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


