from classesForGaussianImport import Atom, utils

scaleFactorExpToBlenderAtoms = .005

def color255to1(color255):
    r = color255[0] / 255.
    g = color255[1] / 255.
    b = color255[2] / 255.
    a = color255[3] / 255.
    color1 = (r, g, b, a)
    return color1

BLENDER = False

if BLENDER:
    import bpy

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
    'Wasserstoff': color255to1(COLOR_RGB255['Schwarz25']),
    'Sauerstoff': color255to1(COLOR_RGB255['Rot100'])
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
        self.size = utils.dictRadiusExp[utils.getKeyFromValue(name, Atom.symbolToName)]
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
    vertices = []
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
        # self.gewichtung = gewichtung
        # self.vertices = []
        # for point in VERTICES_CUBE:
        #    self.vertices.append([point[0] + self.a, point[1] + self.b, point[2] + self.c])

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


def drawBond(startPoint, endPoint, name='default Bond', collection='Collection'):
    if BLENDER:
        bondData = bpy.data.curves.new(name + '_data', 'CURVE')
        bondObject = bpy.data.objects.new(name + '_object', bondData)
        bond = bondData.splines.new('NURBS')  # [POLY, BEZIER, BSPLINE, CARDINAL, NURBS]
        startPoint = startPoint + [1]
        endPoint = endPoint + [1]
        bond.points.add(1)
        bond.points[0].co = startPoint
        bond.points[1].co = endPoint
        bond.use_endpoint_u = True
        bpy.data.collections[collection].objects.link(bondObject)
        bondData.dimensions = '3D'
        bondObject.data.bevel_depth = 0.02
    else:
        print(name + '_data')
        print(name + '_object')
        startPoint = startPoint + [1]
        endPoint = endPoint + [1]
        print(startPoint)
        print(endPoint)



