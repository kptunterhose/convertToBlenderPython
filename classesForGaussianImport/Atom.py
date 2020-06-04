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
  'Cs': 'Caesium', 'Ba': 'Barium', 'La': 'Lanthan', 'Ce': 'Cer', 'Pr': 'Praseodym', 'Nd': 'Neodym', 'Pm': 'Promethium',
  'Sm': 'Samarium', 'Eu': 'Europium', 'Gd': 'Gadolinium', 'Tb': 'Terbium', 'Dy': 'Dysprosium', 'Ho': 'Holmium',
  'Er': 'Erbium', 'Tm': 'Thulium', 'Yb': 'Ytterbium', 'Lu': 'Lutetium', 'Hf': 'Hafnium', 'Ta': 'Tantal', 'W': 'Wolfram',
  'Re': 'Rhenium', 'Os': 'Osmium', 'Ir': 'Iridium', 'Pt': 'Platin', 'Au': 'Gold', 'Hg': 'Quecksilber', 'Tl': 'Thallium',
  'Pb': 'Blei', 'Bi': 'Bismut', 'Po': 'Polonium', 'At': 'Astat', 'Rn': 'Radon',
  'Fr': 'Francium', 'Ra': 'Radium', 'Ac': 'Acinium', 'Th': 'Thorium', 'Pa': 'Proactinium', 'U': 'Uran',
  'Np': 'Neptunium', 'Pu': 'Plutonium', 'Am': 'Americium', 'Cm': 'Curium', 'Bk': 'Berkelium', 'Cf': 'Californium',
  'Es': 'Einsteinium', 'Fm': 'Fermium', 'Md': 'Mendelevium', 'No': 'Nobelium', 'Lr': 'Lawrencium',
  'Rf': 'Rutherfordium', 'Db': 'Dubnium', 'Sg': 'Seaborgium', 'Bh': 'Bohrium', 'Hs': 'Hassium', 'Mt': 'Meitnerium',
  'Ds': 'Darmstadtium', 'Rg': 'Roentgenium', 'Cn': 'Copernicum', 'Nh': 'Nihonium', 'Fl': 'Flerovium', 'Mc': 'Moscovium',
  'Lv': 'Livermorium', 'Ts': 'Tenness', 'Og': 'Oganesson'
}
symbolToNumber = {
  'H': 1, 'He': 2,
  'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
  'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
  'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
  'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
  'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
  'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
  'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65,
  'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77,
  'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86,
  'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,  'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
  'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108,
  'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118
}


def getNameFromSymbol(symbol):
  name =symbolToName[symbol]
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

