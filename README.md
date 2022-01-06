# convertToBlenderPython

### Description

Main: creates python scripts, which creates 3D molekuels from 
quantum chemistry calculations and visualize it in Blender. 
As Input files this script accepts job files from Gaussian (.com),
output files from irc calculation in Gaussian (.log), output from
maxima calculation from AMOLQC and most XYZ files (.xyz).

It also can create NBO calculation for molpro interpreter, or
create wavefunctions for IboView in PSI4. Orbital data of IboView
can also transfers to Blender Python scripts to creat this orbitals
in Blender over the molecuels.


### Dependencies

>+ default: 
>  + ```os, sys, getopt, numpy, json, queue```
>+ irc path with 2 log files (reaction with transition product):
>  + ```scipy```
>+ energy graph of irc path 
>  + ```matplotlib```
>+ creation of qr codes for trigger images:
>  + ```qrcode```
>+ quantum chemistry calculations
>  + ```psi4```

### usage of QuantumChemistryToBlender.py

 ```
 python QuantumChemistryToBlender.py [options] [input files]
 
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
```

### Examples
```QuantumChemistryToBlender.py -c h2o.com
QuantumChemistryToBlender.py -i Reaktion_irc.log
QuantumChemistryToBlender.py -d -1 -i Irc_teil-1.log Irc_teil-2.log
QuantumChemistryToBlender.py -o 9,11-13 # TODO überprüfen mit dem -
```
### TODO
> + in merged file: 
> -x, -w, -m function with 2 files