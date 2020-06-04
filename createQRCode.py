from qrcode.main import QRCode
import os
# quality
LOW = 1
MEDIUM = 0
QUALITY = 3
HIGH = 2


def makeCode(nameShort, nameLong):
    qr = QRCode(
        version=1,
        error_correction=HIGH,
        box_size=20,
        border=4,
    )
    qr.add_data(nameShort)
    qr.make(fit=True)

    img = qr.make_image(fill_color="black", back_color="white")
    img = img.convert('L')
    #img = img.convert('RGB')
    img.save(nameLong + '.png', 'PNG')

listOfAlkans = {
    'Methan': 'CH4',
    'Ethan': 'H3CCH3',
    'Propan': 'H3CCH2CH3',
    'Butan': 'H3CCH2CH2CH3',
    'Pentan': 'H3CCH2CH2CH2CH3',
    'Hexan': 'H3CCH2CH2CH2CH2CH3',
    'Heptan': 'H3CCH2CH2CH2CH2CH2CH3',
    'Octan': 'H3CCH2CH2CH2CH2CH2CH2CH3',
    'Nonan': 'H3CCH2CH2CH2CH2CH2CH2CH2CH3',
    'Decan': 'H3CCH2CH2CH2CH2CH2CH2CH2CH2CH3',
}
listOfAlkens = {
    'Ethen': 'H2CCH2',
    'Propen': 'H2CCHCH3',
    'But-1-en': 'H2CCHCH2CH3',
    'But-2-en': 'H3CCHCHCH3',
    'Pent-1-en': 'H2CCHCH2CH2CH3',
    'Pent-2-en': 'H3CCHCHCH2CH3',
    'Hex-1-en': 'H2CCHCH2CH2CH2CH3',
    'Hex-2-en': 'H3CCHCHCH2CH2CH3',
    'Hex-3-en': 'H3CCH2CHCHCH2CH3',
    'Hept-1-en': 'H2CCHCH2CH2CH2CH2CH3',
    'Hept-2-en': 'H3CCHCHCH2CH2CH2CH3',
    'Hept-3-en': 'H3CCH2CHCHCH2CH2CH3',
    'Oct-1-en': 'H2CCHCH2CH2CH2CH2CH2CH3',
    'Oct-2-en': 'H3CCHCHCH2CH2CH2CH2CH3',
    'Oct-3-en': 'H3CCH2CHCHCH2CH2CH2CH3',
    'Oct-4-en': 'H3CCH2CH2CHCHCH2CH2CH3',
    'Non-1-en': 'H2CCHCH2CH2CH2CH2CH2CH2CH3',
    'Non-2-en': 'H3CCHCHCH2CH2CH2CH2CH2CH3',
    'Non-3-en': 'H3CCH2CHCHCH2CH2CH2CH2CH3',
    'Non-4-en': 'H3CCH2CH2CHCHCH2CH2CH2CH3',
    'Dec-1-en': 'H2CCHCH2CH2CH2CH2CH2CH2CH2CH3',
    'Dec-2-en': 'H3CCHCHCH2CH2CH2CH2CH2CH2CH3',
    'Dec-3-en': 'H3CCH2CHCHCH2CH2CH2CH2CH2CH3',
    'Dec-4-en': 'H3CCH2CH2CHCHCH2CH2CH2CH2CH3',
    'Dec-5-en': 'H3CCH2CH2CH2CHCHCH2CH2CH2CH3'
}
listOfAlkins = {
    'Ethen': 'HCCH',
    'Propen': 'HCCCH3',
    'But-1-en': 'HCCCH2CH3',
    'But-2-en': 'H3CCCCH3',
    'Pent-1-en': 'HCCCH2CH2CH3',
    'Pent-2-en': 'H3CCCCH2CH3',
    'Hex-1-en': 'HCCCH2CH2CH2CH3',
    'Hex-2-en': 'H3CCCCH2CH2CH3',
    'Hex-3-en': 'H3CCH2CCCH2CH3',
    'Hept-1-en': 'HCCCH2CH2CH2CH2CH3',
    'Hept-2-en': 'H3CCCCH2CH2CH2CH3',
    'Hept-3-en': 'H3CCH2CCCH2CH2CH3',
    'Oct-1-en': 'HCCCH2CH2CH2CH2CH2CH3',
    'Oct-2-en': 'H3CCCCH2CH2CH2CH2CH3',
    'Oct-3-en': 'H3CCH2CCCH2CH2CH2CH3',
    'Oct-4-en': 'H3CCH2CH2CCCH2CH2CH3',
    'Non-1-en': 'HCCCH2CH2CH2CH2CH2CH2CH3',
    'Non-2-en': 'H3CCCCH2CH2CH2CH2CH2CH3',
    'Non-3-en': 'H3CCH2CCCH2CH2CH2CH2CH3',
    'Non-4-en': 'H3CCH2CH2CCCH2CH2CH2CH3',
    'Dec-1-en': 'HCCCH2CH2CH2CH2CH2CH2CH2CH3',
    'Dec-2-en': 'H3CCCCH2CH2CH2CH2CH2CH2CH3',
    'Dec-3-en': 'H3CCH2CCCH2CH2CH2CH2CH2CH3',
    'Dec-4-en': 'H3CCH2CH2CCCH2CH2CH2CH2CH3',
    'Dec-5-en': 'H3CCH2CH2CH2CCCH2CH2CH2CH3'
}
listOfAlkohols = {
    'Methanol': 'H3COH',
    'Ethanol': 'H3CCH2OH',
    'Propan-1-ol': 'H3CCH2CH2OH',
    'Propan-2-ol': 'H3CCHOHCH3',
    'Butan-1-ol': 'H3CCH2CH2CH2OH',
    'Butan-2-ol': 'H3CCHOHCH2CH3',
    'Pentan-1-ol': 'H3CCH2CH2CH2CH2OH',
    'Pentan-2-ol': 'H3CCHOHCH2CH2CH3',
    'Pentan-3-ol': 'H3CCH2CHOHCH2CH3',
    'Hexan-1-ol': 'H3CCH2CH2CH2CH2CH2OH',
    'Hexan-2-ol': 'H3CCHOHCH2CH2CH2CH3',
    'Hexan-3-ol': 'H3CCH2CHOHCH2CH2CH3',
    'Heptan-1-ol': 'H3CCH2CH2CH2CH2CH2CH2OH',
    'Heptan-2-ol': 'H3CCHOHCH2CH2CH2CH2CH3',
    'Heptan-3-ol': 'H3CCH2CHOHCH2CH2CH2CH3',
    'Heptan-4-ol': 'H3CCH2CH2CHOHCH2CH2CH3',
    'Octan-1-ol': 'H3CCH2CH2CH2CH2CH2CH2CH2OH',
    'Octan-2-ol': 'H3CCHOHCH2CH2CH2CH2CH2CH3',
    'Octan-3-ol': 'H3CCH2CHOHCH2CH2CH2CH2CH3',
    'Octan-4-ol': 'H3CCH2CH2CHOHCH2CH2CH2CH3',
    'Nonan-1-ol': 'H3CCH2CH2CH2CH2CH2CH2CH2CH2OH',
    'Nonan-2-ol': 'H3CCHOHCH2CH2CH2CH2CH2CH2CH3',
    'Nonan-3-ol': 'H3CCH2CHOHCH2CH2CH2CH2CH2CH3',
    'Nonan-4-ol': 'H3CCH2CH2CHOHCH2CH2CH2CH2CH3',
    'Nonan-5-ol': 'H3CCH2CH2CH2CHOHCH2CH2CH2CH3',
    'Decan-1-ol': 'H3CCH2CH2CH2CH2CH2CH2CH2CH2CH2OH',
    'Decan-2-ol': 'H3CCHOHCH2CH2CH2CH2CH2CH2CH2CH3',
    'Decan-3-ol': 'H3CCH2CHOHCH2CH2CH2CH2CH2CH2CH3',
    'Decan-4-ol': 'H3CCH2CH2CHOHCH2CH2CH2CH2CH2CH3',
    'Decan-5-ol': 'H3CCH2CH2CH2CHOHCH2CH2CH2CH2CH3'
}
listOfAldehyde = {
    'Methanal': 'OCH2',
    'Ethanal': 'OCHCH3',
    'Propanal': 'OCHCH2CH3',
    'Butanal': 'OCHCH2CH2CH3',
    'Pentanal': 'OCHCH2CH2CH2CH3',
    'Hexanal': 'OCHCH2CH2CH2CH2CH3',
    'Heptanal': 'OCHCH2CH2CH2CH2CH2CH3',
    'Octanal': 'OCHCH2CH2CH2CH2CH2CH2CH3',
    'Nonanal': 'OCHCH2CH2CH2CH2CH2CH2CH2CH3',
    'Decanal': 'OCHCH2CH2CH2CH2CH2CH2CH2CH2CH3'
}
listOfKetons = {
    'Propan-2-on': 'H3CCOCH3',
    'Butan-2-on': 'H3CCOCH2CH3',
    'Pentan-2-on': 'H3CCOCH2CH2CH3',
    'Pentan-3-on': 'H3CCH2COCH2CH3',
    'Hexan-2-on': 'H3CCOCH2CH2CH2CH3',
    'Hexan-3-on': 'H3CCH2COCH2CH2CH3',
    'Heptan-2-on': 'H3CCOCH2CH2CH2CH2CH3',
    'Heptan-3-on': 'H3CCH2COCH2CH2CH2CH3',
    'Heptan-4-on': 'H3CCH2CH2COCH2CH2CH3',
    'Octan-2-on': 'H3CCOCH2CH2CH2CH2CH2CH3',
    'Octan-3-on': 'H3CCH2COCH2CH2CH2CH2CH3',
    'Octan-4-on': 'H3CCH2CH2COCH2CH2CH2CH3',
    'Nonan-2-on': 'H3CCOCH2CH2CH2CH2CH2CH2CH3',
    'Nonan-3-on': 'H3CCH2COCH2CH2CH2CH2CH2CH3',
    'Nonan-4-on': 'H3CCH2CH2COCH2CH2CH2CH2CH3',
    'Nonan-5-on': 'H3CCH2CH2CH2COCH2CH2CH2CH3',
    'Decan-2-on': 'H3CCOCH2CH2CH2CH2CH2CH2CH2CH3',
    'Decan-3-on': 'H3CCH2COCH2CH2CH2CH2CH2CH2CH3',
    'Decan-4-on': 'H3CCH2CH2COCH2CH2CH2CH2CH2CH3',
    'Decan-5-on': 'H3CCH2CH2CH2COCH2CH2CH2CH2CH3',
}
listOfCarbonsaeuren = {
    'Methansäure': 'HCOOH',
    'Ethansäure': 'H3CCOOH',
    'Propansäure': 'H3CCH2COOH',
    'Butansäure': 'H3CCH2CH2COOH',
    'Pentansäure': 'H3CCH2CH2CH2COOH',
    'Hexansäure': 'H3CCH2CH2CH2CH2COOH',
    'Heptansäure': 'H3CCH2CH2CH2CH2CH2COOH',
    'Octansäure': 'H3CCH2CH2CH2CH2CH2CH2COOH',
    'Nonansäure': 'H3CCH2CH2CH2CH2CH2CH2CH2COOH',
    'Decansäure': 'H3CCH2CH2CH2CH2CH2CH2CH2CH2COOH'
}
listOfEthers = {
    'Dimethylether': 'H3COCH3',
    'Ethylmethylether': 'H3CCH2OCH3',
    'Diethylether': 'H3CCH2OCH2CH3',
    'Propylmethylether': 'H3CCH2CH2OCH3',
    'Propylethylether': 'H3CCH2CH2OCH2CH3',
    'Dipropylether': 'H3CCH2CH2OCH2CH2CH3',
    'Butylmethylether': 'H3CCH2CH2CH2OCH3',
    'Butylethylether': 'H3CCH2CH2CH2OCH2CH3',
    'Butylpropylether': 'H3CCH2CH2CH2OCH2CH2CH3',
    'Dibutylether': 'H3CCH2CH2CH2OCH2CH2CH2CH3'
}
listOfEsters = {
    'Methansäuremethylester': 'HCOOCH3',
    'Methansäureethylester': 'HCOOCH2CH3',
    'Methansäurepropylester': 'HCOOCH2CH2CH3',
    'Ethansäuremethylester': 'H3CCOOCH3',
    'Ethensäureethylester': 'H3CCOOCH2CH3',
    'Ethansäurepropylester': 'H3CCOOCH2CH3CH3',
    'Propansäuremethylester': 'H3CCH2COOCH3',
    'Propansäureetyhlester': 'H3CCH2COOCH2CH3',
    'Propansäurepropylester': 'H3CCH2COOCH2CH2CH3',
}

dictOfAll = {
    'Alkane': listOfAlkans,
    'Alkene': listOfAlkens,
    'Alkine': listOfAlkins,
    'Alkohole': listOfAlkohols,
    'Aldehyde': listOfAldehyde,
    'Ketone': listOfKetons,
    'Carbonsäuren': listOfCarbonsaeuren,
    'Ether': listOfEthers,
    'Ester': listOfEsters
}

for klasse in dictOfAll:
    if not os.path.exists(klasse):
        os.makedirs(klasse)
    for lN in dictOfAll[klasse]:
        saveTo = klasse + '/' + lN
        makeCode(dictOfAll[klasse][lN], saveTo)




