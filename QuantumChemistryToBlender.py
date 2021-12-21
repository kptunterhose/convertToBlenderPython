#!/usr/bin/python3

import os
import sys
import getopt
import numpy as np
import json
import queue as q
PSI4LOAD = False
MATPLOTLIPLOAD = False
SCIPYLOAD = False
QRCODELOAD = False
VERBOSE = False
IRC_ONEFILE = True
IRC_DIRECTION = 1  # 1 or -1
SORTandRENAME = False
CSV = 'CSV'
LOG = 'LOG'
MOLDEN = 'MOLDEN'
MOLPRO = 'MOLPRO'
MP4 = 'MP4'
XYZ = 'XYZ'
WEBM = 'WEBM'
FILE_TYPES = {
    CSV: 'csv',
    LOG: 'log',
    MOLPRO: 'com',
    MOLDEN: 'molden',
    MP4: 'mp4',
    XYZ: 'xyz',
    WEBM: 'webm'
}
MODULENOTFOUND_CONDA_BASE = 'Das Modul "{}" kann nicht gefunden werden. Schnelle Hilfe wenn es benötigt wird: > conda activate base'
MODULENOTFOUND_CONDA_P4DEV = 'Das Modul "{}" kann nicht gefunden werden. Schnelle Hilfe wenn es benötigt wird: > conda activate p4dev'


def importStuff():
    global PSI4LOAD, MATPLOTLIPLOAD, SCIPYLOAD, QRCODELOAD
    global MODULENOTFOUND_CONDA_BASE, MODULENOTFOUND_CONDA_P4DEV
    global QRCode, plt, animation, minimize, psi4
    try:
        from qrcode.main import QRCode
        QRCODELOAD = True
    except ModuleNotFoundError:
        print(MODULENOTFOUND_CONDA_BASE.format('qrcode'))

    try:
        import matplotlib
        from matplotlib import pyplot as plt
        from matplotlib import animation
        matplotlib.use('TkAgg')
        MATPLOTLIPLOAD = True
    except ModuleNotFoundError:
        print(MODULENOTFOUND_CONDA_BASE.format('matplotlib'))

    try:
        from scipy.optimize import minimize
        SCIPYLOAD = True
    except ModuleNotFoundError:
        print(MODULENOTFOUND_CONDA_BASE.format('scipy'))

    try:
        import psi4
        PSI4LOAD = True
    except ModuleNotFoundError:
        print(MODULENOTFOUND_CONDA_P4DEV.format('psi4'))
        if VERBOSE:
            print('''
            Could not load/find PSI 4, so i can't do any quantum chemistry jobs.
            If you already installd psi 4, you have to load it:
            > conda activate p4dev
            And rerun this script again. If you don't need PSI 4 everything 
            is fine.
            Wenn du IRC Pfade mit mehr als 1 .log Datei in Blender Skripte
            umwandeln willst muss wieder das base Modul geladen werden:
            > conda activate base
            
            
            Mit:
            > conda info --envs
            können die Umgebungen aufgelistet werden 
            base -> scipy, matplotlib, qrcode
            p4dev -> psi4
            ''')


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


def printHelp():
    print("""
    Gebauch von QuantumChemistryToBlender:
    (python) QuantumChemistryToBlender.py [Optionen] [Input-Dateien]
    Bei PSI4 jobs ist es notwendig python vorne weg zu schreiben, da es ansonsten /usr/bin/python3 nutzt
    und nicht die der conda Umgebung (conda info --envs) 
    Python Pfad: witch python
    
    Beispiele:
    QuantumChemistryToBlender.py -c h2o.com 
    QuantumChemistryToBlender.py -i Reaktion_irc.log
    QuantumChemistryToBlender.py -d -1 -i Irc_teil-1.log Irc_teil-2.log
    QuantumChemistryToBlender.py -o 9,11-13 # TODO überprüfen
    
    
    -h        --help            zeigt diese Hilfe an
    -v        --verbose         Zeigt mehr Output im Terminal an
    -c        --com             wandelt .com Datei (Gaussian Job File) in ein Python Skript für Blender um
    -i        --irc             wandelt .log Datei einer IRC Pfad Rechnung in ein Python Skript für Blender um
    -d [1/-1] --direction       Option bei -i mit zwei .log Dateien, gibt die Richtung der Reaktion aus der zweiten .log Datei an: 1 für vorwärts, -1 für rückwärts
    -g        --graph           Erstellt einen animierten Graph: Energie vs. Reaktionskoordinate 
    -m        --molpro          wandelt .log Datei einer IRC Pfad Rechnung in Molpro Job Files (.com) um
    -x        --xyz             wandelt .log Datei einer IRC Pfad Rechnung in xyz Geometrie Dateien für PSI4 Rechnungen um
    -w        --wellenfunktion  erstellt aus .log Datei einer IRC Pfad Rechnung Wellenfunktionen, die in einer Molden Datei gespeichert werden. Diese können mit IboView geöffent und bearbeitet werden.
    -a        --max             wandelt .max/.out Datei einer Maximarechnung aus AMOLQC in ein Python Skript für Blender um
    -o [n]    --orbitals        erstellt aus .orb Dateien aus IboView-Mod Python Skripte für Blender um die Orbitale darzustellen
    """)


BLENDER_CONTENT_HEADER = '''
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
BLENDER_CONTENT_HEADER_ORBITALS = '''
import bpy

CollectionName = 'Orbital'

dOrbs = '''
BLENDER_CONTENT_COM = '''
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
    'Kohlenstoff': color255to1(COLOR_RGB255['Schwarz100']),
    'Wasserstoff': color255to1(COLOR_RGB255['Schwarz10']),
    'Sauerstoff': color255to1(COLOR_RGB255['Rot100']),
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
BLENDER_CONTENT_PATHDATA = '''
"""
pathData = """
'''
BLENDER_CONTENT_IRC = '''
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
BLENDER_CONTENT_MAX = '''
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
BLENDER_CONTENT_ORBITALS = '''
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
'''

MOLPRO_HEADER = '''
memory,200,m;

geometry={
'''
MOLPRO_JOB_DETAILS = '''
}

basis=tzvpp
{df-rks,pbe}
{ibba,iboexp=2, bonds=1}
{put,xml,%s.xml; keepspherical; nosort; skipvirt}
'''


def makeIRCGraph(fileName):
    #fileName: csv mit # Energie Reaktionskoordinate als spalten
    if MATPLOTLIPLOAD:
        outFileName = fileName.split('.')[0]
        extension = fileName.split('.')[-1]

        data = list()
        if extension == FILE_TYPES[CSV]:
            with open(fileName, 'r') as f:
                flines = f.readlines()

            for line in flines:
                line = line.split(';')
                data.append((float(line[1].replace(',', '.')), float(line[2].replace(',', '.'))))
        elif extension == FILE_TYPES[LOG]:
            inFile = open(fileName, 'r')
            while 'Summary of reaction path following' not in inFile.readline():
                pass
            inFile.readline()  # dashes
            inFile.readline()  # Titles of Axes
            line = inFile.readline()  # first line of data
            while '---' not in line:
                line = line.split()
                data.append((float(line[1]), float(line[2])))
                line = inFile.readline()

        xData = list()
        yData = list()
        for dat in data:
            yData.append(dat[0])
            xData.append(dat[1])

        x = xData
        y = yData

        fig, ax = plt.subplots()
        fig.patch.set_alpha(0.)
        fig.set_figwidth(10)
        fig.set_figheight(5)
        plt.grid(True)
        plt.xticks([-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10], fontsize=15)
        plt.yticks(fontsize=15)
        ax.set_xlabel('Reaktionskoordinate', fontsize=15)
        ax.set_ylabel('Energie', fontsize=15)

        # line, = ax.plot(x, y, color='r', linewidth=3)

        point, = ax.plot([], [], marker='+', markersize=7000, markeredgewidth=2, color='r', linewidth=10)

        plt.plot(x, y, linestyle='-', color='b',
                 linewidth=2)  # linestyle: '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'

        def update(num, x, y, line):
            line.set_data(x[:num], y[:num])
            line.axes.axis([min(xData) * 1.1, max(xData) * 1.1, min(yData) * 1.1, max(yData) - min(yData) * .1])
            return line,

        def updateDot(num, x, y, line):
            point.set_data(x[num], y[num])
            point.axes.axis([min(xData) * 1.1, max(xData) * 1.1, min(yData) * 1.1, max(yData) - min(yData) * .1])
            return point,

        ani = animation.FuncAnimation(fig, updateDot, len(x), fargs=[x, y, line],
                                      interval=50, blit=True)

        ani.save(outFileName + '.' + FILE_TYPES[MP4],
                 codec='png',
                 dpi=300,
                 fps=24,
                 savefig_kwargs={'transparent': True})
        os.system('ffmpeg -i {mp4} -c:v libvpx -pix_fmt yuva420p -auto-alt-ref 0 -b:v 4M -y {webm}'.format(mp4=outFileName+'.'+FILE_TYPES[MP4], webm=outFileName+'.'+FILE_TYPES[WEBM]))
        plt.show()
    pass


def getBesucht(adjazenzlist, start):
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


class MaxFile(object):

    def __init__(self, fileName):
        self.inFileName = fileName
        self.outFileName = fileName.split('.')[0] + '_blender.py'

    def inFileToString(self, nMax=1, subNMax=1):
        with open(self.inFileName, 'r') as inFile:
            line = inFile.readline()
            self.nNuclei = int(line.split()[0])
            self.linesToWrite = list()
            for i in range(self.nNuclei):
                words = inFile.readline().split()
                line = ''
                for word in words:
                    line += word + ' '
                line += '\n'
                self.linesToWrite.append(line)
            self.nMaxFound = int(inFile.readline().split()[0])
            maximumData = self.getMaxData(inFile)
            electronLines = maximumData[nMax][subNMax]['lines']
        for line in electronLines:
            self.linesToWrite.append(line())

    def getMaxData(self, openFile, nMax=1):
        data = openFile.readline()
        maxData = dict()
        for i in range(nMax):
            maxData[i+1] = dict()
        while data != '':
            data = data.split()
            actualMax = int(data[1])
            actualSubMax = int(data[2])
            valuesPsiSq = float(data[4])
            nFound = int(data[6])
            nElectrons = int(openFile.readline())
            listOfLines = list()
            for e in range(nElectrons):
                coords = openFile.readline()
                lineToWrite = 'e ' + coords
                listOfLines.append(lineToWrite)
            maxData[actualMax][actualSubMax] = {
                'valuesPsiSquared': valuesPsiSq,
                'nMaxFound': nFound,
                'lines': listOfLines
            }
            data = openFile.readline()
        return maxData

    def makePythonOutfile(self):
        self.inFileToString()
        with open(self.outFileName, 'w') as outFile:
            outFile.write(BLENDER_CONTENT_HEADER)
            for line in self.linesToWrite:
                outFile.write(line)
            outFile.write(BLENDER_CONTENT_MAX)
        print('Speichere Blender Python Skript unter: ' + self.outFileName)


class ComFile(object):

    def __init__(self, fileName):
        self.optionsLines = list()
        self.jobLine = str()
        self.nameLine = str()
        self.multiLine = str()
        self.geometryLines = list()
        self.offSet = list()
        if not fileName:
            fileName = input('please enter file name of gaussian file:\n')
        self.fileName = fileName
        self.outFileName = fileName.split('.')[0] + '_blender.py'

    def readComFile(self):
        countLinesAfterJob = 0
        for i in range(len(self.rawLines)):
            if self.rawLines[i][0] == '%':
                self.optionsLines.append(self.rawLines[i])
            elif self.rawLines[i][0] == '#':
                self.jobLine = self.rawLines[i]
            else:
                countLinesAfterJob += 1
                if countLinesAfterJob == 2:
                    self.nameLine = self.rawLines[i]
                elif countLinesAfterJob == 4:
                    self.multiLine = self.rawLines[i]
                elif countLinesAfterJob > 4:
                    self.geometryLines.append(self.rawLines[i])

    def moveToCenter(self):
        self.calcCenterOfAtoms()
        self.moveAllAtoms()

    def calcCenterOfAtoms(self):
        sum = [0, 0, 0]
        for line in self.geometryLines:
            line = line.split()
            if line:
                sum[0] += float(line[1])
                sum[1] += float(line[2])
                sum[2] += float(line[3])
        for i in range(len(sum)):
            sum[i] = sum[i] / len(lines)
        self.offSet = sum

    def moveAllAtoms(self):
        newLines = []
        for i in range(len(self.geometryLines)):
            line = ''
            words = self.geometryLines[i].split()
            if words:
                line += words[0] + ' '
                line += str(float(words[1]) - offset[0]) + ' '
                line += str(float(words[2]) - offset[1]) + ' '
                line += str(float(words[3]) - offset[2]) + '\n'
                newLines.append(line)
        self.geometryLines = newLines

    def getListOfFiles(self, fileType='com'):
        fileNames = []
        temp = os.listdir()
        for i in range(len(temp)):
            name = temp[i].split('.')
            if name[-1] == fileType:
                fileNames.append(temp[i])
        self.listOfFiles = fileNames

    def makePythonOutfile(self):
        with open(self.fileName, 'r') as inFile:
            self.rawLines = inFile.readlines()
        self.readComFile()
        with open(self.outFileName, 'w') as outFile:
            outFile.write(BLENDER_CONTENT_HEADER)
            for line in self.geometryLines:
                outFile.write(line)
            outFile.write(BLENDER_CONTENT_COM)
        print('Speichere Blender Pyhton Skript unter: ' + self.outFileName)


class IrcFile(object):

    def __init__(self, listOfFiles):
        self.coords4RC = dict()
        self.reactionCoordinates = list()
        self.startConfig = list()  # müsste list sein
        self.linesGeometry = list()
        self.outFileList = list()
        if len(listOfFiles) == 1:
            inFileNamePartOne = listOfFiles[0]
            inFileNamePartTwo = str()
        elif len(listOfFiles) == 2:
            IRC_ONEFILE = False
            inFileNamePartOne = listOfFiles[0]
            inFileNamePartTwo = listOfFiles[1]
        else:
            inFileNamePartOne = input('Bitte geben den Pfad der .log Datei der irc Pfad Rechnung an:\n')
        self.inFileNamePartOne = inFileNamePartOne
        self.inFileNamePartTwo = inFileNamePartTwo
        self.baseFileName = inFileNamePartOne.split('.')[0]
        self.outFileName = self.baseFileName + '_blender.py'
        self.readIrcLog()

    def readIrcLog(self):
        with open(self.inFileNamePartOne, 'r') as inFile:
            self.allLines = inFile.readlines()
        self.getStartGeometry()

    def getStartLines(self):
        startPoints = []
        inputOrientation = []
        inputOrientation4Point = {}
        for i in range(len(self.allLines)):
            if 'Point Number' in self.allLines[i]:
                startPoints.append(i)
            elif 'Input orientation:' in self.allLines[i]:
                inputOrientation.append(i)
        for p in startPoints:
            maxinputPoint = 0
            for inputO in inputOrientation:
                if (inputO < p and inputO > maxinputPoint):
                    maxinputPoint = inputO
            inputOrientation4Point[p] = maxinputPoint
        return startPoints, inputOrientation4Point

    def getStartGeometry(self):
        startP, dict4input = self.getStartLines()

        for point in dict4input:
            pointNR, pathNr, reactionCoordinate = self.getInfos(point)
            self.coords4RC[reactionCoordinate] = self.getKoords(dict4input[point])
            self.reactionCoordinates.append(reactionCoordinate)
        self.reactionCoordinates.sort()
        self.startConfig = self.coords4RC[self.reactionCoordinates[0]]

    def getInfos(self, startLine):
        words = self.allLines[startLine].split()
        pointNR = int(words[2])
        pathNR = int(words[5])
        if pointNR == 0:
            reactionCoordinate = 0
        else:
            reactionCoordinate = float(self.allLines[startLine + 2].split()[-1])
        if pathNR == 1:
            reactionCoordinate *= IRC_DIRECTION
        elif pathNR == 2:
            reactionCoordinate *= -IRC_DIRECTION
        return pointNR, pathNR, reactionCoordinate

    def getKoords(self, startLine):
        if VERBOSE:
            print('reading koords with startline: ' + str(startLine))
        koords = dict()
        lineNumber = startLine + 5
        while '---' not in self.allLines[lineNumber]:
            words = self.allLines[lineNumber].split()
            centerNumber = int(words[0])
            atomicNumber = int(words[1])
            x = float(words[3])
            y = float(words[4])
            z = float(words[5])
            koords[centerNumber] = [atomicNumber, [x, y, z]]
            lineNumber += 1
        return koords

    def makeStartGeometryLines(self):
        for i in self.startConfig:
            lineGeo = numberToSymbol[self.startConfig[i][0]]
            lineGeo += '   '
            lineGeo += str(self.startConfig[i][1][0]) + '  '
            lineGeo += str(self.startConfig[i][1][1]) + '  '
            lineGeo += str(self.startConfig[i][1][2]) + '\n'
            self.linesGeometry.append(lineGeo)

    def getReaktionsPfadLines(self):
        linesOutput = list()
        for j in self.reactionCoordinates:
            lineTemp = 'reactions Coordinate: '
            lineTemp += str(j)
            lineTemp += ' Energie: \n'
            linesOutput.append(lineTemp)
            for centerNumber in self.coords4RC[j]:
                lineTemp = str(centerNumber)
                lineTemp += ' '
                lineTemp += str(self.coords4RC[j][centerNumber][1][0])
                lineTemp += ' '
                lineTemp += str(self.coords4RC[j][centerNumber][1][1])
                lineTemp += ' '
                lineTemp += str(self.coords4RC[j][centerNumber][1][2])
                lineTemp += '\n'
                linesOutput.append(lineTemp)
        return linesOutput

    def getReaktionsPfadLinesFurReaktionsKoordinate(self, reactionsKoordinate):
        linesOutput = list()
        for j in range(len(self.coords4RC[reactionsKoordinate])):
            jj = j + 1
            lineTemp = str(numberToSymbol[self.startConfig[jj][0]])
            lineTemp += ' '
            lineTemp += str(self.coords4RC[reactionsKoordinate][jj][1][0])
            lineTemp += ' '
            lineTemp += str(self.coords4RC[reactionsKoordinate][jj][1][1])
            lineTemp += ' '
            lineTemp += str(self.coords4RC[reactionsKoordinate][jj][1][2])
            lineTemp += '\n'
            linesOutput.append(lineTemp)
        return linesOutput

    def makePythonOutfile(self):
        with open(self.outFileName, 'w') as outFile:
            outFile.write(BLENDER_CONTENT_HEADER)
            self.makeStartGeometryLines()
            for line in self.linesGeometry:
                outFile.write(line)
            outFile.write(BLENDER_CONTENT_PATHDATA)
            for line in self.getReaktionsPfadLines():
                outFile.write(line)
            if not IRC_ONEFILE:
                if SCIPYLOAD:
                    if VERBOSE:
                        print('starte mit zweiter Datei')

                    compareGeometry = self.coords4RC[0]
                    oldReactionsCoordinate = 0
                    for reaktionsKoordinate in self.coords4RC:
                        if reaktionsKoordinate > oldReactionsCoordinate:
                            oldReactionsCoordinate = reaktionsKoordinate
                            compareGeometryPartTwo = self.coords4RC[reaktionsKoordinate]
                    reaktionsKoordinatenShift += abs(reaktionsKoordinate)

                    with open(inFileNamePartTwo, 'r') as inFile:
                        self.allLines = inFile.readlines()
                    getStartGeometry()
                    compareGeometryPartTwo = self.coords4RC[0]
                    oldReactionsCoordinate = 0
                    for reaktionsKoordinate in self.coords4RC:
                        if reaktionsKoordinate < oldReactionsCoordinate:
                            oldReactionsCoordinate = reaktionsKoordinate
                            compareGeometryPartTwo = coords4RCTwo[reaktionsKoordinate]
                    reaktionsKoordinatenShift += abs(reaktionsKoordinate)

                    # calc rot und move matrices --> min distance of coords of fist 6 atoms

                    startPoints = list()
                    targetPoints = list()
                    for i in range(1, 7):
                        startPoints.append(compareGeometryPartTwo[i][1])
                        targetPoints.append(compareGeometry[i][1])
                    rotTrans = RotTrans2Congruent(startPoints=startPoints,
                                                  targetPoints=targetPoints)
                    rot, trans = rotTrans.calc_rot_trans()

                    # make rot und trans
                    # xyz = rotTrans.matrix_mal_vector(rot, xyz)
                    # xyz = rotTrans.translations(trans, xyz)
                    # print(xyz)

                    for line in self.getReaktionsPfadLines():
                        try:
                            # print(len(line.split()))
                            n = int(line.split()[0])
                            x = float(line.split()[1])
                            y = float(line.split()[2])
                            z = float(line.split()[3])
                            vector = [x, y, z]
                            vector = rotTrans.matrix_mal_vektor(rot, vector)
                            vector = rotTrans.translation(trans, vector)
                            line = str(n) + ' ' + str(vector[0]) + ' ' + str(vector[1]) + ' ' + str(vector[2]) + '\n'
                            outFile.write(line)
                        except ValueError:
                            shift = float(line.split()[2]) + reactions_coordinate_shift
                            line = 'reactions Coordinate: '
                            line += str(shift)
                            line += ' Energie: \n'
                            # 'change the reactionscoordinate '
                            outFile.write(line)
                else:
                    print(
                        'konnte scipy nicht laden, nutze base Umgebung (conda activate base), und somit auch nicht die Molekülgeometrien vergleichen')
            outFile.write(BLENDER_CONTENT_IRC)
            print('Speichere Blender Skript unter: ' + self.outFileName)

    def makeGeometryOutfiles(self, format=MOLPRO):
        if os.path.exists(self.baseFileName):
            print('Der Ordner (' + self.baseFileName + ') existiert.')
            overrideInput = input('Wollen Sie den Inhalt überschreiben? [j]/n ')
            if overrideInput.lower() in ('n', 'no', 'nein'):
                print('Skript wird beendet.')
                print('Hinweis: Die Dateien werden in dem Ordner mit dem Namen der .log Datei (ohne Endung) gespeichert')
                exit(1)
        else:
            os.mkdir(self.baseFileName)
            print('Erstelle Output-Ordner: ' + self.baseFileName)
        i = 0

        for reaktionsKoordinate in self.reactionCoordinates:
            reactionsInt = str(i)
            while len(reactionsInt) < 4:
                reactionsInt = '0' + reactionsInt
            i += 1
            outFileName = self.baseFileName + '/' + str(reactionsInt) + '.' + FILE_TYPES[format]
            self.outFileList.append(outFileName)
            with open(outFileName, 'w') as outFile:
                if format == MOLPRO:
                    outFile.write(MOLPRO_HEADER)
                for line in self.getReaktionsPfadLinesFurReaktionsKoordinate(reaktionsKoordinate):
                    outFile.write(line)
                if format == MOLPRO:
                    outFile.write(MOLPRO_JOB_DETAILS % reactionsInt)
            print('Speichere Output-Datei unter: ' + outFileName)

        if SORTandRENAME:
            listOfFiles = os.listdir(self.baseFileName)
            listOfFiles.sort(key=self.sortFunction)
            for i in range(len(listOfFiles)):
                os.rename(self.baseFileName + '/' + listOfFiles[i],
                          self.baseFileName + '/irc-' + str(i) + '.' + FILE_TYPES[format])

    def makeWavefunctionMoldenOutfiles(self):
        global PSI4LOAD
        if PSI4LOAD:
            psi4.set_memory('500 MB')
            psi4.set_options({'PARALLEL': True,
                              'reference': 'rhf'})
            geometries = dict()
            energies = dict()
            wavefunctions = dict()
            for file in self.outFileList:
                print(file.split('.')[-1])
                print(XYZ)
                if file.split('.')[-1] != FILE_TYPES[XYZ]:
                    print('Entferne ' + file + ' von der Liste, da es keine xyz Datei ist.')
                    self.outFileList.remove(file)
                else:
                    with open(file) as f:
                        geometries[file] = psi4.geometry(f.read())
                    print('Berechne Wellenfunktion für ' + file)
                    energies[file], wavefunctions[file] = psi4.energy('mp2/cc-pvdz',
                                                                      molecule=geometries[file],
                                                                      return_wfn=True)
                    outFileName = file.split('.')[0] + '.' + FILE_TYPES[MOLDEN]
                    print('Speichere Molden Datei unter: ' + outFileName)
                    psi4.molden(wavefunctions[file], outFileName)
        else:
            print('kann PSI4 nicht laden')
            print('Versuche Skript mit:')
            print('> python QuantumChemistryToBlender.py -w irc.log')

    def sortFunction(self, value):
        return float(value[:-4])


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


class OrbitalFile(object):

    def __init__(self, onlyList):
        self.dOrbs = dict()
        self.onlyList = onlyList
        self.basePath = './'
        self.outFileName = self.basePath + 'drawOrbitals_blender.py'
        self.listOfOrbitalNames = list()
        listOfNames = os.listdir(self.basePath)
        for name in listOfNames:
            if name.split('.')[-1] == 'orb':
                self.listOfOrbitalNames.append(name)
        self.readOrbitalFromFiles()

    def readOrbitalFromFiles(self):
        for orbitalName in self.listOfOrbitalNames:
            add = False
            if len(self.onlyList) == 0:
                add = True
            else:
                for only in self.onlyList:
                    if 'o' + only in orbitalName:
                        add = True
            if add:
                self.dOrbs[orbitalName] = dict()

                print('Lese Orbital: ' + orbitalName)

                orb = Orbital(self.basePath + orbitalName)
                orb.makePosNegList()
                orb.makePosOrbs()
                orb.makeNegOrbs()
                orb.doFacesStuff()
                dictOfFacesPos, dictOfVertsPos = orb.getPostiveFacesAndVertsDict()
                dictOfFacesNeg, dictOfVertsNeg = orb.getNegativeFacesAndVertsDict()

                if not orb.testStuff():
                    print('Irgendwas mit den Flächen und Vertices passt nicht')

                for nSurface in dictOfFacesPos:
                    self.dOrbs[orbitalName]['pos-' + str(nSurface)] = {
                        'verts': dictOfVertsPos[nSurface],
                        'faces': dictOfFacesPos[nSurface]
                    }
                for nSurface in dictOfFacesNeg:
                    self.dOrbs[orbitalName]['neg-' + str(nSurface)] = {
                        'verts': dictOfVertsNeg[nSurface],
                        'faces': dictOfFacesNeg[nSurface]
                    }

    def makePythonOutput(self):
        with open(self.outFileName, 'w') as outFile:
            outFile.write(BLENDER_CONTENT_HEADER_ORBITALS)
            outFile.write(json.dumps(self.dOrbs))
            outFile.write(BLENDER_CONTENT_ORBITALS)
        print('Speichere Blender Python Skript unter: ' + self.outFileName)


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
            dictOfFaces[surface] = list()
            for face in faces:
                for point in face:
                    if point in dictOfSurfaces[surface]:
                        dictOfFaces[surface].append(list(face))
                        break
            dictOfVerts[surface] = list()
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
        surfaces = dict()
        nSurface = 0
        while listOfKnots:
            s = getBesucht(knoten, listOfKnots[0])
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
    todo = sys.argv
    if len(todo) == 1:
        printHelp()
    try:
        opts, inputFiles = getopt.getopt(sys.argv[1:], 'hvcid:gmxwao:', ['help',
                                                                         'verbose',
                                                                         'com',
                                                                         'irc',
                                                                         'direction=',
                                                                         'graph',
                                                                         'molpro',
                                                                         'xyz',
                                                                         'wellenfunktion',
                                                                         'max',
                                                                         'orbitals='])
    except getopt.GetoptError:
        printHelp()
        sys.exit(2)
    # checking nach globalen optionen
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            printHelp()
        elif opt in ('-v', '--verbose'):
            VERBOSE = True
        elif opt in ('-d', '--direction'):
            if arg == '1':
                IRC_DIRECTION = 1
            else:
                IRC_DIRECTION = -1

    # checking nach jop type
    for opt, arg in opts:
        if opt in ('-c', '--com'):
            for inputFile in inputFiles:
                molekule = ComFile(inputFile)
                molekule.makePythonOutfile()
        elif opt in ('-i', '-m', '-x', '-w', '--irc', '--molpro', 'xyz', 'wellenfunktion'):
            importStuff()
            irc = IrcFile(inputFiles)
            if opt in ('-i', '--irc'):
                irc.makePythonOutfile()
            elif opt in ('-m', '--molpro'):
                irc.makeGeometryOutfiles(format=MOLPRO)
            elif opt in ('-x', '--xyz'):
                irc.makeGeometryOutfiles(format=XYZ)
            elif opt in ('-w', '--wellenfunktion'):
                irc.makeGeometryOutfiles(format=XYZ)
                irc.makeWavefunctionMoldenOutfiles()
        elif opt in ('-g', '--graph'):
            importStuff()
            makeIRCGraph(inputFiles[0])
        elif opt in ('-a', '--max'):
            for f in inputFiles:
                maxFile = MaxFile(fileName=f)
                maxFile.makePythonOutfile()
        elif opt in ('-o', '--orbitals'):
            orbitalFile = OrbitalFile(arg)
            orbitalFile.makePythonOutput()



