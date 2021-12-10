#!/usr/bin/python3

import sys
import os
try:
    import psi4
    PSI4LOAD = True
    CALCWF = False
except ModuleNotFoundError:
    print('''
    Could not load/find PSI 4, so i can't do any quantum chemistry jobs.
    If you already installd psi 4, you have to load it:
    > conda activate p4dev
    And rerun this script again. If you don't need PSI 4 everything 
    is fine.
    ''')
    PSI4LOAD = False



ONEFILE = True
SORTandRENAME = False
xyzOUTPUT = True

symbolToName = {
    'H': 'Wasserstoff', 'He': 'Helium',
    'Li': 'Lithium', 'Be': 'Berillium', 'B': 'Bor', 'C': 'Kohlenstoff', 'N': 'Stickstoff', 'O': 'Sauerstoff',
    'F': 'Flour', 'Ne': 'Neon',
    'Na': 'Natrium', 'Mg': 'Magnesium', 'Al': 'Aluminium', 'Si': 'Silizium', 'P': 'Phosphor', 'S': 'Schwefel',
    'Cl': 'Chlor', 'Ar': 'Argon',
    'K': 'Kalium', 'Ca': 'Kalzium', 'Sc': 'Scandium', 'Ti': 'Titan', 'V': 'Vanadium', 'Cr': 'Chrom', 'Mn': 'Mangan',
    'Fe': 'Eisen', 'Co': 'Kobalt', 'Ni': 'Nickel', 'Cu': 'Kupfer', 'Zn': 'Zink', 'Ga': 'Gallium', 'Ge': 'Germanium',
    'As': 'Arsen', 'Se': 'Selen', 'Br': 'Brom', 'Kr': 'Krypton',
    'Rb': 'Rubidium', 'Sr': 'Strontium', 'Y': 'Yttrium', 'Zr': 'Zirconium', 'Nb': 'Niob', 'Mo': 'Molybdän',
    'Tc': 'Techneticum', 'Ru': 'Ruthenium', 'Rh': 'Rhodium', 'Pd': 'Palladium', 'Ag': 'Silber', 'Cd': 'Cadmium',
    'In': 'Indium', 'Sn': 'Zinn', 'Sb': 'Antimon', 'Te': 'Tellur', 'I': 'Jod', 'Xe': 'Xenon',
    'Cs': 'Caesium', 'Ba': 'Barium', 'La': 'Lanthan', 'Ce': 'Cer', 'Pr': 'Praseodym', 'Nd': 'Neodym',
    'Pm': 'Promethium',
    'Sm': 'Samarium', 'Eu': 'Europium', 'Gd': 'Gadolinium', 'Tb': 'Terbium', 'Dy': 'Dysprosium', 'Ho': 'Holmium',
    'Er': 'Erbium', 'Tm': 'Thulium', 'Yb': 'Ytterbium', 'Lu': 'Lutetium', 'Hf': 'Hafnium', 'Ta': 'Tantal',
    'W': 'Wolfram',
    'Re': 'Rhenium', 'Os': 'Osmium', 'Ir': 'Iridium', 'Pt': 'Platin', 'Au': 'Gold', 'Hg': 'Quecksilber',
    'Tl': 'Thallium',
    'Pb': 'Blei', 'Bi': 'Bismut', 'Po': 'Polonium', 'At': 'Astat', 'Rn': 'Radon',
    'Fr': 'Francium', 'Ra': 'Radium', 'Ac': 'Acinium', 'Th': 'Thorium', 'Pa': 'Proactinium', 'U': 'Uran',
    'Np': 'Neptunium', 'Pu': 'Plutonium', 'Am': 'Americium', 'Cm': 'Curium', 'Bk': 'Berkelium', 'Cf': 'Californium',
    'Es': 'Einsteinium', 'Fm': 'Fermium', 'Md': 'Mendelevium', 'No': 'Nobelium', 'Lr': 'Lawrencium',
    'Rf': 'Rutherfordium', 'Db': 'Dubnium', 'Sg': 'Seaborgium', 'Bh': 'Bohrium', 'Hs': 'Hassium', 'Mt': 'Meitnerium',
    'Ds': 'Darmstadtium', 'Rg': 'Roentgenium', 'Cn': 'Copernicum', 'Nh': 'Nihonium', 'Fl': 'Flerovium',
    'Mc': 'Moscovium',
    'Lv': 'Livermorium', 'Ts': 'Tenness', 'Og': 'Oganesson'
}
symbolToNumber = {
    'H': 1, 'He': 2,
    'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
    'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
    'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47,
    'Cd': 48,
    'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
    'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65,
    'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76,
    'Ir': 77,
    'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
    'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
    'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107,
    'Hs': 108,
    'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}
numberToSymbol = {
    1: 'H', 2: 'He',
    3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
    11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
    19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
    30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
    37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag',
    48: 'Cd',
    49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe',
    55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb',
    66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os',
    77: 'Ir',
    78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
    87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk',
    98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh',
    108: 'Hs',
    109: 'Mt', 110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'
}


def getStartLines(lines):
    startPoints = []
    inputOrientation = []
    inputOrientation4Point = {}
    for i in range(len(lines)):
        if 'Point Number' in lines[i]:
            startPoints.append(i)
        elif 'Input orientation:' in lines[i]:
            inputOrientation.append(i)
    for p in startPoints:
        maxinputPoint = 0
        for inputO in inputOrientation:
            if (inputO < p and inputO > maxinputPoint):
                maxinputPoint = inputO
        inputOrientation4Point[p] = maxinputPoint
    return startPoints, inputOrientation4Point


def getStartGeometry(allLines):
    startP, dict4input = getStartLines(allLines)

    for point in dict4input:
        pointNR, pathNr, reactionCoordinate = getInfos(point)
        coords4RC[reactionCoordinate] = getKoords(dict4input[point])

    reactionCoordinates = []
    for p in coords4RC:
        reactionCoordinates.append(p)
    reactionCoordinates.sort()

    startConfig = coords4RC[reactionCoordinates[0]]
    return startConfig, reactionCoordinates


def getInfos(startLine):
    words = allLines[startLine].split()
    pointNR = int(words[2])
    pathNR = int(words[5])
    if pointNR == 0:
        reactionCoordinate = 0
    else:
        reactionCoordinate = float(allLines[startLine + 2].split()[-1])
    if pathNR == 1:
        reactionCoordinate *= DIRECTION
    elif pathNR == 2:
        reactionCoordinate *= -DIRECTION
    return pointNR, pathNR, reactionCoordinate


def getKoords(startLine):
    if VERBOSE:
        print('reading koords with startline: ' + str(startLine))
    koords = {}
    lineNumber = startLine + 5
    while '---' not in allLines[lineNumber]:
        words = allLines[lineNumber].split()
        centerNumber = int(words[0])
        atomicNumber = int(words[1])
        x = float(words[3])
        y = float(words[4])
        z = float(words[5])
        koords[centerNumber] = [atomicNumber, [x, y, z]]
        lineNumber += 1
    return koords


def getStartGeometryLines(startConfig):
    linesGeometry = []
    for i in startConfig:
        lineGeo = numberToSymbol[startConfig[i][0]]
        lineGeo += '   '
        lineGeo += str(startConfig[i][1][0]) + '  '
        lineGeo += str(startConfig[i][1][1]) + '  '
        lineGeo += str(startConfig[i][1][2]) + '\n'
        linesGeometry.append(lineGeo)
    return linesGeometry


def getReaktionsPfadLines(reactionsCoord, coords4ReactionsPath, startConfig):
    linesOutput = []
    for j in range(len(coords4ReactionsPath[reactionsCoord])):
        jj = j+1
        lineTemp = str(numberToSymbol[startConfig[jj][0]])
        lineTemp += ' '
        lineTemp += str(coords4ReactionsPath[reactionsCoord][jj][1][0])
        lineTemp += ' '
        lineTemp += str(coords4ReactionsPath[reactionsCoord][jj][1][1])
        lineTemp += ' '
        lineTemp += str(coords4ReactionsPath[reactionsCoord][jj][1][2])
        lineTemp += '\n'
        linesOutput.append(lineTemp)
    return linesOutput


def sortFunction(value):
    return float(value[:-4])


content1 = '''
memory,200,m;

geometry={
'''

content2 = '''
}

basis=tzvpp
{df-rks,pbe}
{ibba,iboexp=2, bonds=1}
{put,xml,%s.xml; keepspherical; nosort; skipvirt}
'''


if __name__ == '__main__':
    if '-wf' in sys.argv:
        CALCWF = True
        sys.argv.remove('-wf')
    ### noch eine inputabfrage wie wf in sys.argv CALCWF = True
    if len(sys.argv) == 2:
        inFileName_part_one = sys.argv[1]
    elif len(sys.argv) == 3:
        ONEFILE = False
        inFileName_part_one = sys.argv[1]
        inFileName_part_two = sys.argv[2]
    else:
        inFileName_part_one = 'Veresterung_irc.log'
        #inFileName_part_one = input('please enter file name of gaussian file:\n')

    VERBOSE = False
    DIRECTION = 1  # 1 OR -1
    coords4RC = {}

    baseFileName = inFileName_part_one.split('.')[0]

    inFile = open(inFileName_part_one, 'r')
    allLines = inFile.readlines()
    startConfig, rcs = getStartGeometry(allLines=allLines)

    if os.path.exists(baseFileName):
        print('folder already exists')
        print('program exec now')
        override_input = input('want to override y/n?')
        if override_input.lower() == 'n':
            exit(1)
    else:
        os.mkdir(baseFileName)
        print('create folder: ' + baseFileName)

    ### quick an dirty
    i = 0
    outFileList = []
    for reactionsCoord in rcs:
        reactionsInt = str(i)
        while len(reactionsInt) < 4:
            reactionsInt = '0' + reactionsInt
        i += 1
        outFileName = baseFileName + '/' + str(reactionsInt) + '.'
        if xyzOUTPUT:
            outFileName += 'xyz'
        else:
            outFileName += 'com'
        outFileList.append(outFileName)
        outFile = open(outFileName, 'w')
        if not xyzOUTPUT:
            outFile.write(content1)
        for line in getReaktionsPfadLines(reactionsCoord, coords4RC, startConfig):
            outFile.write(line)
        if not xyzOUTPUT:
            outFile.write(content2 % reactionsInt)
        outFile.close()
        print('write output to ' + outFileName)

    if SORTandRENAME:
        listOfFiles = os.listdir(baseFileName)
        listOfFiles.sort(key=sortFunction)
        for i in range(len(listOfFiles)):
            os.rename(baseFileName + '/' + listOfFiles[i], baseFileName + '/irc-' + str(i) + '.com')

    if PSI4LOAD and xyzOUTPUT and CALCWF:
        # set psi4 stuff
        psi4.set_memory('500 MB')
        psi4.set_options({'PARALLEL': True,
                          'reference': 'rhf'})
        geometries = {}
        energies = {}
        wavefunctions = {}

        for file in outFileList:
            if file.split('.')[-1] != 'xyz':
                print('remove ' + file + ' from list, because it is no xyz file')
                outFileList.remove(file)
            else:
                with open(file) as f:
                    geometries[file] = psi4.geometry(f.read())
                print('calc wfn for ' + file)
                energies[file], wavefunctions[file] = psi4.energy('mp2/cc-pvdz',
                                                                  molecule=geometries[file],
                                                                  return_wfn=True)
                outFileName = file.split('.')[0] + '.molden'
                print('save molden file: ' + outFileName)
                psi4.molden(wavefunctions[file], outFileName)



