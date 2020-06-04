from classesForGaussianImport import Atom, BlenderUtils, utils

scaleFactorExpToBlenderAtoms = .005
scaleFactorSingleBond = .25  # Depth in Blender
maxDistanceForSingleBond = 150.  # Anström


def calcBonds(listOfAtoms):
    bonds = []
    for i in range(len(listOfAtoms)):
        for ii in range(i + 1, len(listOfAtoms)):
            dist = utils.getDistance(listOfAtoms[i].getPosition(), listOfAtoms[ii].getPosition())
            if (dist < maxDistanceForSingleBond / 100.):
                # print(maxDistanceForSingleBond/100.)
                # print(dist)
                bonds.append((listOfAtoms[i], listOfAtoms[ii]))
    return bonds


def readComFile(nameInFile):
    atoms = []
    infile = open(nameInFile, 'r')
    line = infile.readline()
    while line[0] == '%':
        line = infile.readline()
    while line[0] == '#':
        line = infile.readline()
        line = infile.readline()
        line = infile.readline()
        line = infile.readline()
    lines = infile.readlines()

    for i in range(len(lines)):
        words = lines[i].split()
        if len(words) == 4:
            atoms.append(Atom.Atom(float(words[1]), float(words[2]), float(words[3]), words[0]))
    infile.close()
    return atoms


def readFromString(input):
    atoms = []
    lines = input.split('\n')
    for line in lines:
        words = line.split()
        if len(words) == 4:
            atoms.append(Atom.Atom(float(words[1]), float(words[2]), float(words[3]), words[0]))
    return atoms

    pass


def getElements(listOfAtoms):
    elements = set()
    for atom in listOfAtoms:
        elements.add(atom.getSymbol())
    return elements


if __name__ == '__main__':
    atome = readComFile('MeOH_opt.com')
    #atome = readFromString(inputtext)
    elemente = getElements(atome)
    bindungen = calcBonds(atome)
    materialElemente = {}
    atome4Blender = []
    for element in elemente:
        print(element)
        materialElemente[element] = BlenderUtils.Element4Blender(Atom.symbolToName[element])
    atome4Blender = []
    for atom in atome:
        a, b, c = atom.getPosition()
        atome4Blender.append(BlenderUtils.Atom4Blender(a, b, c, materialElemente[atom.getSymbol()]))
    for a in atome4Blender:
        a.draw()
        #a.addMaterial()
        #a.makeSphere()
    for bindung in bindungen:
        name = 'bond_' + str(bindung[0].getSymbol()) + '-' + str(bindung[1].getSymbol())
        BlenderUtils.drawBond(bindung[0].getPosition(), bindung[1].getPosition(), name)

