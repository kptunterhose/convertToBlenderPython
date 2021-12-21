#!/usr/bin/env python3
'''
TODO:
Bindungen wie in ircToBlender
if sys.args[1] == '*'
convert all .com-Files inside the folder
'''
import sys
import os

content1 = '''
import numpy as np
BLENDER = True
if BLENDER:
    import bpy
    
scaleFactorExpToBlenderAtoms = .015
scaleFactorSingleBond = .07  # Depth in Blender
maxDistanceForSingleBond = 150.  # Anström    

startGeo = """
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
    'Lv': 'Livermorium', 'Ts': 'Tenness', 'Og': 'Oganesson',
    'e': 'Elektron'
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
    'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118,
    'e': 0
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

    def __init__(self, x, y, z, identifier, charge=0):
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

    def isElectron(self):
        if self.symbol == 'e':
            return True
        else:
            return False    


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
    'Nh': 0, 'Fl': 0, 'Mc': 0, 'Lv': 0, 'Ts': 0, 'Og': 0,
    'e': 5
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
    'Nh': 0, 'Fl': 0, 'Mc': 0, 'Lv': 0, 'Ts': 0, 'Og': 0,
    'e': 5
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
    'Kohlenstoff': color255to1(COLOR_RGB255['Schwarz100']),
    'Wasserstoff': color255to1(COLOR_RGB255['Schwarz10']),
    'Sauerstoff': color255to1(COLOR_RGB255['Rot100']),
    'Chlor': color255to1(COLOR_RGB255['Grun100']),
    'Bindung': color255to1(COLOR_RGB255['Schwarz50']),
    'Elektron': color255to1(COLOR_RGB255['Petrol100'])
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
    #vertices = []
    element = ''
    label = ''
    gewichtung = 0

    def __init__(self, a, b, c, element):
        self.a = a
        self.b = b
        self.c = c
        self.position = [a, b, c]
        self.element = element
        self.label = self.element.getName() + '_' + str(NATOMS[self.element.getName()])
        NATOMS['Gesamt'] += 1
        NATOMS[self.element.getName()] += 1

    def draw(self, collectionName='Collection', level=5):
        if BLENDER:
            self.mesh = bpy.data.meshes.new('mesh_' + self.label)
            self.mesh.from_pydata(VERTICES_CUBE, [], FACES_CUBE)
            self.obj = bpy.data.objects.new(self.label, self.mesh)
            self.obj.location = self.position
            self.obj.scale = [self.element.getSize(), self.element.getSize(), self.element.getSize()]
            bpy.data.collections[collectionName].objects.link(self.obj)
            self.makeSphere(level=level)
            self.addMaterial()
        else:
            print('male atom')

    def makeSphere(self, level=5):
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


class Bond4Blender(object):
    startPoint = []
    #startPointLine = []
    endPoint = []
    #endPointLine = []
    farbe = COLORATOM['Bindung']
    bondInactiv = False
    name = ''
    collection = 'Collection'


    def __init__(self, startPoint, endPoint, name, collection='Collection'):
        self.startPoint = []
        self.endPoint = []
        self.setStartPoint(startPoint[0], startPoint[1], startPoint[2])
        self.setEndPoint(endPoint[0], endPoint[1], endPoint[2])
        self.name = name
        self.collection = collection
        if BLENDER:
            self.material = bpy.data.materials.new(name=name)
            self.material.use_nodes = True
            self.material.diffuse_color = self.farbe
        else:
            self.material = 'Material'
            print('material erstellt')
            print('Farbe = ' + str(self.farbe))
            

    def createBond(self):
        if BLENDER:
            self.bondData = bpy.data.curves.new(self.name + '_data', 'CURVE')
            self.bondObject = bpy.data.objects.new(self.name + '_object', self.bondData)
            self.bond = self.bondData.splines.new('NURBS')  # [POLY, BEZIER, BSPLINE, CARDINAL, NURBS]
            self.bond.points.add(1)
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
            self.bond.points[0].co = [0, 0, 0, 1]
            self.bond.points[1].co = [0, 0, 1, 1]
            self.bond.use_endpoint_u = True
            bpy.data.collections[self.collection].objects.link(self.bondObject)
            self.bondData.dimensions = '3D'
            self.bondObject.data.bevel_depth = 0.05
            self.bondObject.data.materials.append(self.material)
            self.moveBond()
        else:
            print('draw bonds')

    def moveBond(self):
        phi = np.pi - self.getWinkelZuZAxis(self.startPoint, self.endPoint)
        drehVektor = self.getNormalenVektor(self.startPoint, self.endPoint)
        axisAngle = [phi] + drehVektor
        self.isBondActiv()
        if self.bondInactiv:
            self.bondObject.scale = [0, 0, 0]
        else:
            self.bondObject.scale = self.calcBondScale(self.getBetrag(self.getVerbindungsVektor(self.startPoint, self.endPoint)))
        self.bondObject.location = self.startPoint
        self.bondObject.rotation_mode = 'AXIS_ANGLE'
        self.bondObject.rotation_axis_angle = axisAngle
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
            if listOfAtoms[i].isElectron() or listOfAtoms[ii].isElectron():
                pass
            else:
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


def readFromString(inputGeo):
    atoms = []
    lines = inputGeo.split('\\n')
    centerNumber = 0
    for line in lines:
        if line != '':
            words = line.split()
            if len(words) == 4:
                atoms.append(Atom(float(words[1]), float(words[2]), float(words[3]), words[0]))
    return atoms

def getElements(listOfAtoms):
    elements = set()
    for atom in listOfAtoms:
        elements.add(atom.getSymbol())
    return elements


if __name__ == '__main__':
    # atome = readComFile('MeOH_opt.com')
    atome = readFromString(startGeo)
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
        atome4Blender.append(Atom4Blender(a, b, c, materialElemente[atom.getSymbol()]))
    for a in atome4Blender:
        a.draw()
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
                                               name)))

    for bond in bindungen4Blender:
        bond.drawBond()
'''


def convertComToString(inFile, moveToCenter=True):
    line = inFile.readline()
    while line[0] == '%':
        line = inFile.readline()
    while line[0] == '#':
        line = inFile.readline()
        line = inFile.readline()
        line = inFile.readline()
        line = inFile.readline()
    lines = inFile.readlines()
    inFile.close()
    if moveToCenter:
        offSet = calcCenterOfAtoms(lines)
        lines = moveAllAtoms(lines, offSet)
    return lines


def convertMaxToString(inFile, nMax=1, subNMax=1):
    line = inFile.readline()
    nNuclei = int(line.split()[0])
    linesToWrite = []
    for i in range(nNuclei):
        line = inFile.readline()
        words = line.split()[1:]
        line = ''
        for word in words:
            line += word + '  '
        line += '\n'
        linesToWrite.append(line)
    nMaxFound = int(inFile.readline().split()[0])
    maximumData= getMaxData(inFile, nMaxFound)
    electronLines = maximumData[nMax][subNMax]['lines']
    for line in electronLines:
        linesToWrite.append(line)
    return linesToWrite


def getMaxData(inFile, nMax=1):
    data = inFile.readline()
    maxData = {}
    for i in range(nMax):
        maxData[i+1] = dict()
    while data != '':
        data = data.split()
        acutalMax = int(data[1])
        actualSubMax = int(data[2])
        valuePsiSq = float(data[4])
        nFound = int(data[6])
        nElectrons = int(inFile.readline())
        linesToWrite = []
        for e in range(nElectrons):
            coords = inFile.readline()
            lineToWrite = 'e ' + coords
            linesToWrite.append(lineToWrite)
        maxData[acutalMax][actualSubMax] = {
            'valuePsiSquared': valuePsiSq,
            'nMaxFound': nFound,
            'lines': linesToWrite}
        data = inFile.readline()
    return maxData


def calcCenterOfAtoms(lines):
    sum = [0, 0, 0]
    for line in lines:
        line = line.split()
        if line:
            sum[0] += float(line[1])
            sum[1] += float(line[2])
            sum[2] += float(line[3])
    for i in range(len(sum)):
        sum[i] = sum[i] / len(lines)
    return sum


def moveAllAtoms(lines, offset):
    newLines = []
    for i in range(len(lines)):
        line = ''
        words = lines[i].split()
        if words:
            line += words[0] + ' '
            line += str(float(words[1]) - offset[0]) + ' '
            line += str(float(words[2]) - offset[1]) + ' '
            line += str(float(words[3]) - offset[2]) + '\n'
            newLines.append(line)
    return newLines


def getListOfFiles():
    fileNames = []
    temp = os.listdir()
    for i in range(len(temp)):
        name = temp[i].split('.')
        if name[-1] == 'ref':
            fileNames.append(temp[i])
    return fileNames


def makePythonOutfile(inFileName):
    inFile = open(inFileName, 'r')
    lines = convertMaxToString(inFile)
    inFile.close()
    outFileName = inFileName.split('.')[0] + '_blender.py'
    outFile = open(outFileName, 'w')
    outFile.write(content1)
    for line in lines:
        outFile.write(line)
    outFile.write(content3)
    outFile.close()
    print('write output to ' + outFileName)
    return 0


if __name__ == '__main__':
    if len(sys.argv) == 2:
        inFileName = str(sys.argv[1])
    else:
        inFileName = input('please enter file name of gaussian file:\n')

    if inFileName == '*':
        listOfFiles = getListOfFiles()
        for fileName in listOfFiles:
            makePythonOutfile(fileName)
    else:
        makePythonOutfile(inFileName)

