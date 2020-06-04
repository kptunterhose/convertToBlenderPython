import numpy as np

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
  'La': 226, 'Ce': 210, 'Pr': 247, 'Nd': 206, 'Pm': 205, 'Sm': 238, 'Eu': 231, 'Gd': 233, 'Tb': 225, 'Dy': 228, 'Ho': 226, 'Er': 226, 'Tm': 222, 'Yb': 222, 'Lu': 217,
  'Hf': 208, 'Ta': 200, 'W': 193, 'Re': 188, 'Os': 185, 'Ir': 180, 'Pt': 177, 'Au': 174, 'Hg': 171,
  'Tl': 156, 'Pb': 154, 'Bi': 143, 'Po': 135, 'At': 127, 'Rn': 120,
  'Fr': 0, 'Ra': 0,
  'Ac': 0, 'Th': 0, 'Pa': 0, 'U': 0, 'Np': 0, 'Pu': 0, 'Am': 0, 'Cm': 0, 'Bk':0, 'Cf': 0, 'Es': 0, 'Fm': 0, 'Md': 0, 'No': 0, 'Lr': 0,
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
  'La': 195, 'Ce': 185, 'Pr': 185, 'Nd': 185, 'Pm': 185, 'Sm': 185, 'Eu': 185, 'Gd': 180, 'Tb': 175, 'Dy': 175, 'Ho': 175, 'Er': 175, 'Tm': 175, 'Yb': 175, 'Lu': 175,
  'Hf': 155, 'Ta': 145, 'W': 135, 'Re': 135, 'Os': 130, 'Ir': 135, 'Pt': 135, 'Au': 135, 'Hg': 150,
  'Tl': 190, 'Pb': 180, 'Bi': 160, 'Po': 190, 'At': 0, 'Rn': 0,
  'Fr': 0, 'Ra': 215,
  'Ac': 195, 'Th': 180, 'Pa': 180, 'U': 175, 'Np': 175, 'Pu': 175, 'Am': 175, 'Cm': 0, 'Bk': 0, 'Cf': 0, 'Es': 0, 'Fm': 0, 'Md': 0, 'No': 0, 'Lr': 0,
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

