#!/usr/bin/python3
import sys
import numpy as np
from scipy.optimize import minimize

ONEFILE = True

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


def getReaktionsPfadLines(reactionsCoords, coords4ReactionsPath):
    linesOutput = []
    for j in reactionsCoords:
        lineTemp = 'reactions Coordinate: '
        lineTemp += str(j)
        lineTemp += ' Energie: \n'
        linesOutput.append(lineTemp)
        for centerNumber in coords4ReactionsPath[j]:
            lineTemp = str(centerNumber)
            lineTemp += ' '
            lineTemp += str(coords4ReactionsPath[j][centerNumber][1][0])
            lineTemp += ' '
            lineTemp += str(coords4ReactionsPath[j][centerNumber][1][1])
            lineTemp += ' '
            lineTemp += str(coords4ReactionsPath[j][centerNumber][1][2])
            lineTemp += '\n'
            linesOutput.append(lineTemp)
    return linesOutput


class RotTrans2Congruent(object):
    alt = []
    neu = []
    neu_u = []
    print_output = False
    rot_matrix = [[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]]

    def __init__(self, start_points, target_points, print_output=False):
        self.alt = start_points
        self.neu = target_points
        self.neu_u = target_points
        self.print_output = print_output

    def matrix_multiplikation_3d(self, A, B):
        # dimensionen checken
        C = 3 * [None]
        for i in range(3):
            C[i] = 3 * [0]
            for j in range(3):
                for k in range(3):
                    C[i][k] += A[i][j] * B[j][k]
        return C

    def matrix_mal_vektor(self, M, v):
        n = 3 * [0]
        for i in range(3):
            for j in range(3):
                n[i] += M[i][j] * v[j]
        return n

    def translation(self, transl, point):
        p = 3 * [None]
        if len(transl) == len(point):
            for i in range(len(transl)):
                p[i] = point[i] + transl[i]
        return p

    def make_rot_z(self, alpha):
        co = float(np.cos(alpha))
        si = float(np.sin(alpha))
        rot = [[co, (-1 * si), 0],
               [si, co, 0],
               [0, 0, 1]]
        return rot

    def make_rot_y(self, alpha):
        co = float(np.cos(alpha))
        si = float(np.sin(alpha))
        rot = [[co, 0, si],
               [0, 1, 0],
               [(-1 * si), 0, co]]
        return rot

    def make_rot_x(self, alpha):
        co = float(np.cos(alpha))
        si = float(np.sin(alpha))
        rot = [[1, 0, 0],
               [0, co, (-1 * si)],
               [0, si, co]]
        return rot

    def rot_um_z(self, alpha, point):
        rot_z = self.make_rot_z(alpha)
        n = self.matrix_mal_vektor(rot_z, point)
        return n

    def rot_um_y(self, alpha, point):
        rot_y = self.make_rot_y(alpha)
        n = self.matrix_mal_vektor(rot_y, point)
        return n

    def rot_um_x(self, alpha, point):
        rot_x = self.make_rot_x(alpha)
        n = self.matrix_mal_vektor(rot_x, point)
        return n

    def equation_to_solve_z(self, alpha):
        alt_rot = len(self.alt) * [None]
        for i in range(len(self.alt)):
            alt_rot[i] = self.rot_um_z(alpha, self.alt[i])
        diff = 0
        for i in range(len(self.alt)):
            for j in range(3):
                diff += (self.neu[i][j] - alt_rot[i][j]) ** 2
        return diff

    def equation_to_solve_y(self, alpha):
        alt_rot = len(self.alt) * [None]
        for i in range(len(self.alt)):
            alt_rot[i] = self.rot_um_y(alpha, self.alt[i])
        diff = 0
        for i in range(len(self.alt)):
            for j in range(3):
                diff += (self.neu[i][j] - alt_rot[i][j]) ** 2
        return diff

    def equation_to_solve_x(self, alpha):
        alt_rot = len(self.alt) * [None]
        for i in range(len(self.alt)):
            alt_rot[i] = self.rot_um_x(alpha, self.alt[i])
        diff = 0
        for i in range(len(self.alt)):
            for j in range(3):
                diff += (self.neu[i][j] - alt_rot[i][j]) ** 2
        return diff

    def calc_rot_trans(self, n_steps_max=10):
        # translation of old --> origin
        trans_alt = 3 * [None]
        for i in range(len(self.alt[0])):
            trans_alt[i] = -1 * self.alt[0][i]

        # tralnslation of target --> origin --> back
        trans_neu = 3 * [None]
        trans_neu_revers = 3 * [None]
        for i in range(len(self.neu[0])):
            trans_neu[i] = -1 * self.neu[0][i]
            trans_neu_revers[i] = self.neu[0][i]
        # do translation of both melekules to origin
        for i in range(len(self.alt)):
            self.alt[i] = self.translation(trans_alt, self.alt[i])

        for i in range(len(self.neu)):
            self.neu[i] = self.translation(trans_neu, self.neu[i])

        # init params of omptimization
        alpha_z = []
        alpha_y = []
        alpha_x = []

        # rotation wit max numer of Steps
        for n in range(n_steps_max):
            solution_z = minimize(self.equation_to_solve_y, 0, method='Nelder-Mead')
            alpha_z.append(float(solution_z['x']))
            self.rot_matrix = self.matrix_multiplikation_3d(self.make_rot_z(float(solution_z['x'])), self.rot_matrix)
            for i in range(len(self.alt)):
                self.alt[i] = self.rot_um_z(alpha_z[n], self.alt[i])
            solution_y = minimize(self.equation_to_solve_y, 0, method='Nelder-Mead')
            alpha_y.append(float(solution_y['x']))
            self.rot_matrix = self.matrix_multiplikation_3d(self.make_rot_y(float(solution_y['x'])), self.rot_matrix)
            for i in range(len(self.alt)):
                self.alt[i] = self.rot_um_y(alpha_y[n], self.alt[i])
            solution_x = minimize(self.equation_to_solve_x, 0, method='Nelder-Mead')
            alpha_x.append(float(solution_x['x']))
            self.rot_matrix = self.matrix_multiplikation_3d(self.make_rot_x(float(solution_x['x'])), self.rot_matrix)
            for i in range(len(self.alt)):
                self.alt[i] = self.rot_um_x(alpha_x[n], self.alt[i])

        final_translation = self.translation(trans_neu_revers, self.matrix_mal_vektor(self.rot_matrix, trans_alt))

        # print out results
        if self.print_output:
            print('\nfirst rotation with: ')
            for i in range(len(self.rot_matrix)):
                print('%9.6f %9.6f %9.6f' % (self.rot_matrix[i][0], self.rot_matrix[i][1], self.rot_matrix[i][2]))
            print('\nthen, translation with:')
            for i in range(len(final_translation)):
                print('%9.6f' % final_translation[i])
            print('\n--- finished ---')
        return self.rot_matrix, final_translation


content1 = '''
import numpy as np
BLENDER = True
if BLENDER:
    import bpy

scaleFactorExpToBlenderAtoms = .015
scaleFactorSingleBond = .07  # Depth in Blender
maxDistanceForSingleBond = 160.  # Anström
smoothnessOfAtoms = 3


startGeo = """
'''

content2 = '''
"""
pathData = """
'''

content3 = '''
"""
# Atom.py
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


def getNameFromSymbol(symbol):
    name = symbolToName[symbol]
    return name


class Atom(object):
    number = 0
    symbol = ''
    element = ''
    x = 0
    y = 0
    z = 0
    charge = 0

    def __init__(self, x, y, z, identifier, centerNumber, pathCoordinates, charge=0):
        self.centerNumber = centerNumber
        self.pathCoordinates = pathCoordinates
        self.number = 0
        self.symbol = ''
        self.element = ''
        self.x = 0
        self.y = 0
        self.z = 0
        self.charge = 0
        if isinstance(identifier, int):
            if (identifier <= 118 & identifier > 0):
                self.number = identifier
                self.symbol = self.getSymbolFromNumber()
                self.element = symbolToName[self.symbol]
            else:
                pass
        elif isinstance(identifier, str):
            self.symbol = identifier
            self.number = symbolToNumber[identifier]
            self.element = symbolToName[identifier]
        else:
            pass
        self.x = x
        self.y = y
        self.z = z
        self.charge = charge

    def getSymbolFromNumber(self, number):
        key_list = list(symbolToNumber.keys())
        val_list = list(symbolToNumber.values())
        symbol = key_list[val_list.index(number)]
        return symbol

    def toString(self):
        return (self.symbol + ' at x = ' + str(self.x) + ', y = ' + str(self.y) + ', z = ' + str(self.z))

    def getPosition(self):
        position = [self.x, self.y, self.z]
        return position

    def getSymbol(self):
        return self.symbol

    def reduceReactionsPath(self):
        reducedPath = {}
        for i in self.pathCoordinates:
            reducedPath[i] = self.pathCoordinates[i][self.centerNumber]
        self.pathCoordinates = reducedPath
        pass

    def getReaktionsPath(self):
        self.reduceReactionsPath()
        return self.pathCoordinates





# utils.py
# theoretisch berechnet: https://doi.org/10.1063/1.1712084
dictRadiusTheo = {
    'H': 53, 'He': 31,
    'Li': 167, 'Be': 112, 'B': 87, 'C': 67, 'N': 56, 'O': 48, 'F': 42, 'Ne': 38,
    'Na': 190, 'Mg': 145, 'Al': 118, 'Si': 111, 'P': 98, 'S': 88, 'Cl': 79, 'Ar': 71,
    'K': 243, 'Ca': 194,
    'Sc': 184, 'Ti': 176, 'V': 171, 'Cr': 166, 'Mn': 161, 'Fe': 156, 'Co': 152, 'Ni': 149, 'Cu': 145, 'Zn': 142,
    'Ga': 136, 'Ge': 125, 'As': 114, 'Se': 103, 'Br': 94, 'Kr': 88,
    'Rb': 265, 'Sr': 219,
    'Y': 212, 'Zr': 206, 'Nb': 198, 'Mo': 190, 'Tc': 183, 'Ru': 178, 'Rh': 173, 'Pd': 169, 'Ag': 165, 'Cd': 161,
    'In': 156, 'Sn': 145, 'Sb': 133, 'Te': 123, 'I': 115, 'Xe': 108,
    'Cs': 298, 'Ba': 253,
    'La': 226, 'Ce': 210, 'Pr': 247, 'Nd': 206, 'Pm': 205, 'Sm': 238, 'Eu': 231, 'Gd': 233, 'Tb': 225, 'Dy': 228,
    'Ho': 226, 'Er': 226, 'Tm': 222, 'Yb': 222, 'Lu': 217,
    'Hf': 208, 'Ta': 200, 'W': 193, 'Re': 188, 'Os': 185, 'Ir': 180, 'Pt': 177, 'Au': 174, 'Hg': 171,
    'Tl': 156, 'Pb': 154, 'Bi': 143, 'Po': 135, 'At': 127, 'Rn': 120,
    'Fr': 0, 'Ra': 0,
    'Ac': 0, 'Th': 0, 'Pa': 0, 'U': 0, 'Np': 0, 'Pu': 0, 'Am': 0, 'Cm': 0, 'Bk': 0, 'Cf': 0, 'Es': 0, 'Fm': 0, 'Md': 0,
    'No': 0, 'Lr': 0,
    'Rf': 0, 'Db': 0, 'Sg': 0, 'Bh': 0, 'Hs': 0, 'Mt': 0, 'Ds': 0, 'Rg': 0, 'Cn': 0,
    'Nh': 0, 'Fl': 0, 'Mc': 0, 'Lv': 0, 'Ts': 0, 'Og': 0
}
# experimentell bestimmt https://doi.org/10.1063/1.1725697
dictRadiusExp = {
    'H': 25, 'He': 0,
    'Li': 145, 'Be': 105, 'B': 85, 'C': 70, 'N': 65, 'O': 60, 'F': 50, 'Ne': 0,
    'Na': 180, 'Mg': 150, 'Al': 125, 'Si': 110, 'P': 100, 'S': 100, 'Cl': 100, 'Ar': 0,
    'K': 220, 'Ca': 180,
    'Sc': 160, 'Ti': 140, 'V': 135, 'Cr': 140, 'Mn': 140, 'Fe': 140, 'Co': 135, 'Ni': 135, 'Cu': 135, 'Zn': 135,
    'Ga': 130, 'Ge': 125, 'As': 115, 'Se': 115, 'Br': 115, 'Kr': 0,
    'Rb': 235, 'Sr': 200,
    'Y': 180, 'Zr': 155, 'Nb': 145, 'Mo': 145, 'Tc': 135, 'Ru': 130, 'Rh': 135, 'Pd': 140, 'Ag': 160, 'Cd': 155,
    'In': 155, 'Sn': 145, 'Sb': 145, 'Te': 140, 'I': 140, 'Xe': 0,
    'Cs': 260, 'Ba': 215,
    'La': 195, 'Ce': 185, 'Pr': 185, 'Nd': 185, 'Pm': 185, 'Sm': 185, 'Eu': 185, 'Gd': 180, 'Tb': 175, 'Dy': 175,
    'Ho': 175, 'Er': 175, 'Tm': 175, 'Yb': 175, 'Lu': 175,
    'Hf': 155, 'Ta': 145, 'W': 135, 'Re': 135, 'Os': 130, 'Ir': 135, 'Pt': 135, 'Au': 135, 'Hg': 150,
    'Tl': 190, 'Pb': 180, 'Bi': 160, 'Po': 190, 'At': 0, 'Rn': 0,
    'Fr': 0, 'Ra': 215,
    'Ac': 195, 'Th': 180, 'Pa': 180, 'U': 175, 'Np': 175, 'Pu': 175, 'Am': 175, 'Cm': 0, 'Bk': 0, 'Cf': 0, 'Es': 0,
    'Fm': 0, 'Md': 0, 'No': 0, 'Lr': 0,
    'Rf': 0, 'Db': 0, 'Sg': 0, 'Bh': 0, 'Hs': 0, 'Mt': 0, 'Ds': 0, 'Rg': 0, 'Cn': 0,
    'Nh': 0, 'Fl': 0, 'Mc': 0, 'Lv': 0, 'Ts': 0, 'Og': 0
}


def getKeyFromValue(Value, dictonary):
  key_list = list(dictonary.keys())
  val_list = list(dictonary.values())
  key = key_list[val_list.index(Value)]
  return key


def getDistance(Point3D1, Point3D2):
    distance = 0
    a = (Point3D1[0] - Point3D2[0]) ** 2
    b = (Point3D1[1] - Point3D2[1]) ** 2
    c = (Point3D1[2] - Point3D2[2]) ** 2
    distance = np.sqrt(a + b + c)
    return distance


# BlenderUtils
def color255to1(color255):
    r = color255[0] / 255.
    g = color255[1] / 255.
    b = color255[2] / 255.
    a = color255[3] / 255.
    color1 = (r, g, b, a)
    return color1


VERTICES_CUBE = [(-.5, -.5, -.5),
                 (-.5, .5, -.5),
                 (.5, .5, -.5),
                 (.5, -.5, -.5),
                 (-.5, -.5, .5),
                 (-.5, .5, .5),
                 (.5, .5, .5),
                 (.5, -.5, .5)]
FACES_CUBE = [(0, 1, 2, 3),
              (7, 6, 5, 4),
              (0, 4, 5, 1),
              (1, 5, 6, 2),
              (2, 6, 7, 3),
              (3, 7, 4, 0)]
COLOR = {
    'blau': (0, .3294, .6235, 1),
    'gruen': (.3412, .6706, .1529, 1),
    'rot': (.8, .0275, .1176, 1),
    'bordeaux': (.6314, .0627, .2078, 1),
}
COLOR_RGB255 = {
    'Blau100': (0, 84, 159, 255),
    'Blau75': (64, 127, 183, 255),
    'Blau50': (142, 186, 229, 255),
    'Blau25': (199, 221, 242, 255),
    'Blau10': (232, 241, 250, 255),
    'Schwarz100': (0, 0, 0, 255),
    'Schwarz75': (100, 101, 103, 255),
    'Schwarz50': (156, 158, 159, 255),
    'Schwarz25': (207, 209, 210, 255),
    'Schwarz10': (236, 237, 237, 255),
    'Magenta100': (227, 0, 102, 255),
    'Magenta75': (233, 96, 136, 255),
    'Magenta50': (241, 158, 177, 255),
    'Magenta25': (249, 210, 218, 255),
    'Magenta10': (253, 238, 240, 255),
    'Geld100': (255, 237, 0, 255),
    'Gelb75': (255, 240, 85, 255),
    'Gelb50': (255, 245, 155, 255),
    'Gelb25': (255, 250, 209, 255),
    'Gelb10': (255, 253, 238, 255),
    'Petrol100': (0, 97, 101, 255),
    'Petrol75': (45, 127, 131, 255),
    'Petrol50': (125, 164, 167, 255),
    'Petrol25': (191, 208, 209, 255),
    'Petrol10': (230, 236, 236, 255),
    'Turkies100': (0, 152, 161, 255),
    'Turkies75': (0, 177, 183, 255),
    'Turkies50': (137, 204, 207, 255),
    'Turkies25': (202, 231, 231, 255),
    'Turkies10': (235, 246, 246, 255),
    'Grun100': (87, 171, 39, 255),
    'Grun75': (141, 192, 96, 255),
    'Grun50': (184, 214, 152, 255),
    'Grun25': (211, 235, 206, 255),
    'Grun10': (242, 247, 236, 255),
    'Maigrun100': (189, 205, 0, 255),
    'Maigrun75': (208, 217, 92, 255),
    'Maigrun50': (224, 230, 154, 255),
    'Maigrun25': (240, 243, 208, 255),
    'Maigrun10': (249, 250, 237, 255),
    'Orange100': (246, 168, 0, 255),
    'Orange75': (250, 190, 80, 255),
    'Orange50': (253, 212, 143, 255),
    'Orange25': (254, 234, 201, 255),
    'Orange10': (255, 247, 234, 255),
    'Rot100': (204, 7, 30, 255),
    'Rot75': (216, 92, 65, 255),
    'Rot50': (230, 150, 121, 255),
    'Rot25': (243, 205, 187, 255),
    'Rot10': (250, 235, 227, 255),
    'Bordeaux100': (161, 16, 53, 255),
    'Bordeaux75': (182, 82, 86, 255),
    'Bordeaux50': (205, 139, 135, 255),
    'Bordeaux25': (229, 197, 192, 255),
    'Bordeaux10': (245, 232, 229, 255),
    'Violett100': (97, 33, 88, 255),
    'Violett75': (131, 78, 117, 255),
    'Violett50': (168, 133, 158, 255),
    'Violett25': (210, 192, 205, 255),
    'Violett10': (237, 229, 234, 255),
    'Lila100': (122, 111, 172, 255),
    'Lila75': (155, 145, 193, 255),
    'Lila50': (188, 181, 215, 255),
    'Lila25': (222, 218, 235, 255),
    'Lila10': (242, 240, 247, 255)
}
COLORATOM = {
    'Wasserstoff': color255to1(COLOR_RGB255['Schwarz10']),
    'Lithium': color255to1(COLOR_RGB255['Lila100']),
    'Kohlenstoff': color255to1(COLOR_RGB255['Schwarz100']),
    'Stickstoff': color255to1(COLOR_RGB255['Blau100']),
    'Sauerstoff': color255to1(COLOR_RGB255['Rot100']),
    'Flour': color255to1(COLOR_RGB255['Grun75']),
    'Chlor': color255to1(COLOR_RGB255['Grun100']),
    'Bindung': color255to1(COLOR_RGB255['Schwarz50'])
}
NATOMS = {
    'Gesamt': 0
}


class Element4Blender(object):
    name = ''
    farbe = (0, 0, 0, 0)
    size = 1

    def __init__(self, name):
        if BLENDER:
            self.material = bpy.data.materials.new(name=name)
        else:
            self.material = 'Material'
            print('material erstellt')
        self.name = name
        self.farbe = COLORATOM[name]
        if BLENDER:
            self.material.diffuse_color = self.farbe
        else:
            print('Farbe = ' + str(self.farbe))
        NATOMS[name] = 0
        self.size = dictRadiusExp[getKeyFromValue(name, symbolToName)]
        self.size *= scaleFactorExpToBlenderAtoms

    def getMaterial(self):
        return self.material

    def getName(self):
        return self.name

    def getSize(self):
        return self.size


class Atom4Blender(object):
    a = 0
    b = 0
    c = 0
    position = []
    reactionPath = []
    #vertices = []
    element = ''
    label = ''
    gewichtung = 0

    def __init__(self, a, b, c, element, reactionPath):
        self.reactionPath = reactionPath
        self.a = a
        self.b = b
        self.c = c
        self.position = [a, b, c]
        self.element = element
        self.label = self.element.getName() + '_' + str(NATOMS[self.element.getName()])
        NATOMS['Gesamt'] += 1
        NATOMS[self.element.getName()] += 1
        # self.gewichtung = gewichtung
        #self.vertices = []
        #for point in VERTICES_CUBE:
        #    self.vertices.append([point[0] + self.a, point[1] + self.b, point[2] + self.c])

    def draw(self, collectionName='Collection', level=smoothnessOfAtoms):
        if BLENDER:
            self.mesh = bpy.data.meshes.new('mesh_' + self.label)
            self.mesh.from_pydata(VERTICES_CUBE, [], FACES_CUBE)
            self.obj = bpy.data.objects.new(self.label, self.mesh)
            self.changePosition(self.position)
            self.obj.scale = [self.element.getSize(), self.element.getSize(), self.element.getSize()]
            bpy.data.collections[collectionName].objects.link(self.obj)
            self.makeSphere(level=level)
            self.addMaterial()
        else:
            print('male atom')

    def makeSphere(self, level=smoothnessOfAtoms):
        if BLENDER:
            self.obj.modifiers.new('subd', type='SUBSURF')
            self.obj.modifiers['subd'].levels = level
            self.obj.modifiers['subd'].render_levels = level
        else:
            print('erstelle Kugel mit level = ' + str(level))

    def addMaterial(self):
        if BLENDER:
            self.obj.data.materials.append(self.element.getMaterial())
        else:
            print('Material hinzugefuegt')

    def getLabel(self):
        return self.label

    def getKoordinaten(self):
        return self.a, self.b, self.c

    def getGewichtung(self):
        return self.gewichtung

    def getGewichteteKoordinaten(self):
        gewichteteKoordinaten = [self.gewichtung, self.a, self.b, self.c]
        return gewichteteKoordinaten

    def changePosition(self, pos):
        if BLENDER:
            self.obj.location = pos
        else:
            print('change pos to: ' + str(pos))
        pass

    def animateReactionPath(self):
        offset = 0.
        for rc in self.reactionPath:
            rc = float(rc)
            if rc > offset:
                offset = rc
        for rc in self.reactionPath:
            if float(rc) == 0:
                pass
            else:
                self.changePosition(self.reactionPath[rc])
                frameNumber = offset - float(rc)
                frameNumber *= 10
                if BLENDER:
                    self.obj.keyframe_insert(data_path='location', frame=frameNumber)
                else:
                    print('insert keyframe')
        pass


class Bond4Blender(object):
    startPoint = []
    #startPointLine = []
    endPoint = []
    #endPointLine = []
    farbe = COLORATOM['Bindung']
    bondInactiv = False
    name = ''
    collection = 'Collection'


    def __init__(self, startPoint, endPoint, name, reactionPath, centerNumber1, centerNumber2, collection='Collection'):
        self.startPoint = []
        self.endPoint = []
        self.setStartPoint(startPoint[0], startPoint[1], startPoint[2])
        self.setEndPoint(endPoint[0], endPoint[1], endPoint[2])
        self.name = name
        self.collection = collection
        self.reactionPath = readPathFromString(reactionPath)
        self.centerNumber1 = centerNumber1
        self.centerNumber2 = centerNumber2
        self.isBondActiv()
        if BLENDER:
            self.material = bpy.data.materials.new(name=name)
            self.material.use_nodes = True
        else:
            self.material = 'Material'
            print('material erstellt')
        if BLENDER:
            self.material.diffuse_color = self.farbe
        else:
            print('Farbe = ' + str(self.farbe))

    def createBond(self):
        if BLENDER:
            self.bondData = bpy.data.curves.new(self.name + '_data', 'CURVE')
            self.bondObject = bpy.data.objects.new(self.name + '_object', self.bondData)
            self.bond = self.bondData.splines.new('NURBS')  # [POLY, BEZIER, BSPLINE, CARDINAL, NURBS]
            self.bond.points.add(1)
            #self.bond.points[0].co = self.startPointLine
            #self.bond.points[1].co = self.endPointLine
            self.bond.points[0].co = [0, 0, 0, 1]
            self.bond.points[1].co = [0, 0, 1, 1]
            self.bond.use_endpoint_u = True
            bpy.data.collections[self.collection].objects.link(self.bondObject)
            self.bondData.dimensions = '3D'
            self.bondObject.data.bevel_depth = 0.05
            self.bondObject.data.materials.append(self.material)

    def drawBond(self):
        if BLENDER:
            self.bondData = bpy.data.curves.new(self.name + '_data', 'CURVE')
            self.bondObject = bpy.data.objects.new(self.name + '_object', self.bondData)
            self.bond = self.bondData.splines.new('NURBS')  # [POLY, BEZIER, BSPLINE, CARDINAL, NURBS]
            self.bond.points.add(1)
            #self.bond.points[0].co = self.startPointLine
            #self.bond.points[1].co = self.endPointLine
            self.bond.points[0].co = [0, 0, 0, 1]
            self.bond.points[1].co = [0, 0, 1, 1]
            self.bond.use_endpoint_u = True
            bpy.data.collections[self.collection].objects.link(self.bondObject)
            self.bondData.dimensions = '3D'
            self.bondObject.data.bevel_depth = 0.05
            self.bondObject.data.materials.append(self.material)

            self.moveBond(self.startPoint, self.endPoint)
        else:
            print('draw bonds')

    def moveBond(self, point1, point2):
        self.startPoint = point1
        self.endPoint = point2
        #phiX, phiY = self.getXandYWinkel(point1, point2)
        phi = np.pi - self.getWinkelZuZAxis(point1, point2)
        drehVektor = self.getNormalenVektor(point1, point2)
        axisAngle = [phi] + drehVektor
        self.isBondActiv()
        if self.bondInactiv:
            self.bondObject.scale = [0, 0, 0]
        else:
            self.bondObject.scale = self.calcBondScale(self.getBetrag(self.getVerbindungsVektor(point1, point2)))
        self.bondObject.location = point1
        self.bondObject.rotation_mode = 'AXIS_ANGLE'
        self.bondObject.rotation_axis_angle = axisAngle
        #self.bondObject.rotation_mode = 'XYZ'
        #self.bondObject.rotation_euler = [phiX, phiY, 0]
        #bpy.context.view_layer.objects.active = bpy.data.objects[self.name]
        #bpy.ops.transform.rotate(value=phiX, orienr_axis='X', orient_type='LOCAL')
        #bpy.ops.transform.rotate(value=phiY, orienr_axis='Y', orient_type='LOCAL')
        pass

    def makeAnimation(self):
        offset = 0.
        for rc in self.reactionPath:
            rc = float(rc)
            if rc > offset:
                offset = rc
        for rc in self.reactionPath:
            if float(rc) == 0:
                pass
            else:
                self.moveBond(self.reactionPath[rc][self.centerNumber1], self.reactionPath[rc][self.centerNumber2])
                self.isBondActiv()

                #self.changePosition(self.reactionPath[rc])
                frameNumber = offset - float(rc)
                frameNumber *= 10
                if BLENDER:
                    self.bondObject.keyframe_insert(data_path='location', frame=frameNumber)
                    self.bondObject.keyframe_insert(data_path='rotation_axis_angle', frame=frameNumber)
                    self.bondObject.keyframe_insert(data_path='scale', frame=frameNumber)
                else:
                    print('insert keyframe')
        pass

    def setStartPoint(self, x, y, z):
        self.startPoint = [x, y, z]
        self.startPointLine = [x, y, z, 1]

    # to remove
    def setEndPoint(self, x, y, z):
        self.endPoint = [x, y, z]
        self.endPointLine = [x, y, z, 1]

    # to remove
    def getStartPoint(self):
        return self.startPoint

    # to remove
    def getEndPoint(self):
        return self.endPoint

    def isBondActiv(self):
        dist = getDistance(self.startPoint, self.endPoint)
        if dist > maxDistanceForSingleBond / 100.:
            self.bondInactiv = True
        else:
            self.bondInactiv = False

    def getBetrag(self, vector):
        s = 0
        for v in vector:
            s += v**2
        betrag = np.sqrt(s)
        return betrag

    def getVerbindungsVektor(self, point1, point2):
        vector = []
        if len(point1) == len(point2):
            for i in range(len(point1)):
                vector.append(point1[i] - point2[i])
        return vector

    # solle weg können
    def getXandYWinkel(self, point1, point2):
        phiX = 0
        phiY = 0
        targetVector = self.getVerbindungsVektor(point1, point2)
        try:
            phiX = 1. * (np.pi / 2 - np.arccos(targetVector[1] / self.getBetrag(targetVector)))
            phiY = -1. * (np.pi / 2 - np.arccos(targetVector[0] / self.getBetrag(targetVector)))
        except Exception as e:
            print(e)
        return phiX, phiY

    def getWinkelZuZAxis(self, startVektor, zielVektor):
        phi = 0
        targetVektor = self.getVerbindungsVektor(startVektor, zielVektor)
        try:
            phi = np.arccos(targetVektor[2]/self.getBetrag(targetVektor))
        except Exception as e:
            print(e)
        return phi

    def getNormalenVektor(self, point1, point2):
        vP2 = self.getVerbindungsVektor(point1, point2)
        vP3 = [0, 0, 1]
        normalenVektor = [0, 0, 0]
        normalenVektor[0] = vP2[1] * vP3[2] - vP2[2] * vP3[1]
        normalenVektor[1] = vP2[2] * vP3[0] - vP2[0] * vP3[2]
        normalenVektor[2] = vP2[0] * vP3[1] - vP2[1] * vP3[0]
        betrag = self.getBetrag(normalenVektor)
        for i in range(len(normalenVektor)):
            normalenVektor[i] /= betrag
        return normalenVektor

    def getCoordinatesFromReactionCoordinate(self, reactionCoordinate):
        pass

    def calcBondScale(self, laenge):
        sx = 1/laenge
        sy = 1/laenge
        sz = laenge
        return [sx, sy, sz]

    # überarbeitet, z.t.
    def moveBondAlt(self, startDelta, endDelta):
        # start and endDelat must be like (dx, dy, dz)
        if BLENDER:
            # set object activ
            bpy.context.view_layer.objects.active = bpy.data.objects[self.name]
            # enter edit-mode
            bpy.ops.object.mode_set(mode='EDIT')
            # deselect all vertices
            bpy.ops.curve.select_all(action='DESELECT')
            # select last vertex
            bpy.ops.curve.de_select_last()
            # translate value=(dx, dy, dz)
            bpy.ops.transform.translate(value=endDelta)
            # deselect all vertices
            bpy.ops.curve.select_all(action='DESELECT')
            # select first vertex
            bpy.ops.curve.de_select_first()
            # translate value=(dx, dy, dz)
            bpy.ops.transform.translate(value=startDelta)
            # deselect all vertices
            bpy.ops.curve.select_all(action='DESELECT')
            # get back in object mode
            bpy.ops.object.mode_set(mode='OBJECT')


def getListOfAllBonds(listOfAtoms):
    bonds = {}
    #bonds = []
    for i in range(len(listOfAtoms)):
        for ii in range(i + 1, len(listOfAtoms)):
            bonds[(i+1, ii+1)] = [listOfAtoms[i], listOfAtoms[ii]]
            #bonds.append((listOfAtoms[i], listOfAtoms[ii]))
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
            atoms.append(Atom(float(words[1]), float(words[2]), float(words[3]), words[0]))
    infile.close()
    return atoms


def readFromString(inputGeo, inputPath):
    atoms = []
    lines = inputGeo.split('\\n')
    rc = readPathFromString(inputPath)
    centerNumber = 0
    for line in lines:
        if line != '':
            centerNumber += 1
            words = line.split()
            if len(words) == 4:
                print(centerNumber)
                atoms.append(Atom(float(words[1]), float(words[2]), float(words[3]), words[0], centerNumber, rc))
    return atoms


def readPathFromString(inputPath):
    """Return a dict with path data for all atoms
    key 1 is the reaction coordinate
    key 2 is the center number of the atom
    value is a list of x, y, and z coordinates of the atom
    :param inputPath:
    :return: dict[reactionCoordinate][centerNumber] = [x, y, z]
    """
    # dict im dict:
    # key 1 ist reaktions koordinate
    # key 2 ist die centernumber des atoms
    # value ist [x, y, z]
    lines = inputPath.split('\\n')
    reaktionCordinate = {}
    for line in lines[1:]:
        if line != '':
            words = line.split()
            if words[0] == 'reactions':
                reakCoord = words[2]
                reaktionCordinate[reakCoord] = {}
            else:
                reaktionCordinate[reakCoord][int(words[0])] = [float(words[1]), float(words[2]), float(words[3])]
    return reaktionCordinate


def getElements(listOfAtoms):
    elements = set()
    for atom in listOfAtoms:
        elements.add(atom.getSymbol())
    return elements


if __name__ == '__main__':
    # atome = readComFile('MeOH_opt.com')
    atome = readFromString(startGeo, pathData)
    elemente = getElements(atome)
    bindungen = getListOfAllBonds(atome)
    materialElemente = {}
    atome4Blender = []
    if BLENDER:
        bpy.context.scene.render.engine = 'CYCLES'

    for element in elemente:
        print(element)
        materialElemente[element] = Element4Blender(symbolToName[element])
    atome4Blender = []
    for atom in atome:
        a, b, c = atom.getPosition()
        pfad = atom.getReaktionsPath()
        atome4Blender.append(Atom4Blender(a, b, c, materialElemente[atom.getSymbol()], pfad))
    for a in atome4Blender:
        a.draw()
        a.animateReactionPath()
        #a.addMaterial()
        #a.makeSphere()
    bindungen4Blender = []
    nBindung = 1
    for cNbrs in bindungen:
        name = 'bond_' \
               + str(nBindung) + '_' \
               + str(bindungen[cNbrs][0].getSymbol()) \
               + '-' + str(bindungen[cNbrs][1].getSymbol())
        nBindung += 1
        bindungen4Blender.append((Bond4Blender(bindungen[cNbrs][0].getPosition(),
                                               bindungen[cNbrs][1].getPosition(),
                                               name,
                                               pathData,
                                               cNbrs[0],
                                               cNbrs[1])))

    for bond in bindungen4Blender:
        bond.createBond()
        bond.makeAnimation()
'''

if __name__ == '__main__':
    if len(sys.argv) == 2:
        inFileName_part_one = sys.argv[1]
    elif len(sys.argv) == 3:
        ONEFILE = False
        inFileName_part_one = sys.argv[1]
        inFileName_part_two = sys.argv[2]
    else:
        # inFileName = 'Veresterung_irc.log'
        inFileName_part_one = input('please enter file name of gaussian file:\n')

    VERBOSE = False
    DIRECTION = 1  # 1 OR -1
    coords4RC = {}

    inFile = open(inFileName_part_one, 'r')
    allLines = inFile.readlines()
    startConfig, rcs = getStartGeometry(allLines=allLines)

    outFileName = inFileName_part_one.split('.')[0] + '_blender.py'
    outFile = open(outFileName, 'w')

    outFile.write(content1)
    for line in getStartGeometryLines(startConfig):
        outFile.write(line)

    outFile.write(content2)

    for line in getReaktionsPfadLines(rcs, coords4RC):
        outFile.write(line)

    if not ONEFILE:
        if VERBOSE:
            print('start with second file')

        compare_geometry = coords4RC[0]
        old_reactions_coordinate = 0
        for reactions_coordinate in coords4RC:
            if reactions_coordinate > old_reactions_coordinate:
                old_reactions_coordinate = reactions_coordinate
                compare_geometry = coords4RC[reactions_coordinate]
        reactions_coordinate_shift = abs(reactions_coordinate)
        coords4RC = {}
        in_file_part_two = open(inFileName_part_two, 'r')
        allLines = in_file_part_two.readlines()

        startConfig_part_two, rcs_part_two = getStartGeometry(allLines=allLines)
        compare_geometry_part_two = coords4RC[0]
        old_reactions_coordinate = 0
        for reactions_coordinate in coords4RC:
            if reactions_coordinate < old_reactions_coordinate:
                old_reactions_coordinate = reactions_coordinate
                compare_geometry_part_two = coords4RC[reactions_coordinate]
        reactions_coordinate_shift += abs(reactions_coordinate)
        # calc rot und move matrixes --> min distance of coords of first 6 atoms
        start_points = []
        for i in range(1, 7):
            start_points.append(compare_geometry_part_two[i][1])
        target_points = []
        for i in range(1, 7):
            target_points.append(compare_geometry[i][1])

        test = RotTrans2Congruent(start_points, target_points)
        rot, trans = test.calc_rot_trans()

        # make rot and trans
        # xyz = test.matrix_mal_vector(rot, xyz)
        # xyz = test.translation(trans, xyz)
        # print xyz
        for line in getReaktionsPfadLines(rcs_part_two, coords4RC):
            try:
                #print(len(line.split()))
                n = int(line.split()[0])
                x = float(line.split()[1])
                y = float(line.split()[2])
                z = float(line.split()[3])
                vector = [x, y, z]
                vector = test.matrix_mal_vektor(rot, vector)
                vector = test.translation(trans, vector)
                line = str(n) + ' ' + str(vector[0]) + ' ' + str(vector[1]) + ' ' + str(vector[2]) + '\n'
                outFile.write(line)
            except ValueError:
                shift = float(line.split()[2]) + reactions_coordinate_shift
                line = 'reactions Coordinate: '
                line += str(shift)
                line += ' Energie: \n'
                #'change the reactionscoordinate '
                outFile.write(line)

    outFile.write(content3)
    outFile.close()
    print('write output to ' + outFileName)

'''
listOfAtoms = set()



exit(3)
startP, dict4input = getStartLines(allLines)

for point in dict4input:
    pointNR, pathNr, reactionCoordinate = getInfos(point)
    coords4RC[reactionCoordinate] = getKoords(dict4input[point])

reactionCoordinates = []
for p in coords4RC:
    reactionCoordinates.append(p)
reactionCoordinates.sort()

startConfig = [reactionCoordinates[0], coords4RC[reactionCoordinates[0]]]

for k in coords4RC[reactionCoordinates[0]]:
    print(k, coords4RC[reactionCoordinates[0]][k])
    listOfAtoms.add(k)


'''
