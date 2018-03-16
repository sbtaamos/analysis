# import generic python modules
# hello
import argparse
import operator
from operator import itemgetter
import sys, os, shutil, itertools
import os.path

###################################################################################################

# RETRIEVE USER INPUTS

###################################################################################################

###################################################################################################

# create parser

###################################################################################################

version_nb = "0.1.0"
parser = argparse.ArgumentParser(prog='cluster_prot', usage='', add_help=False, formatter_class=argparse.RawDescriptionHelpFormatter, description=\
'''
**********************************************
v''' + version_nb + '''
author: Sarah-Beth Amos (sarah-beth.amos@bioch.ox.ac.uk)
git: 
**********************************************


	
[ DESCRIPTION ]

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
  
''')

import mdtraj as md
import numpy as np
import mdanalysis as mda
import numpy as np
import MDAnalysis
#from MDAnalysis.core.parallel.distances import distance_array
from MDAnalysis.analysis.distances import distance_array, self_distance_array
import networkx as nx
import subprocess
import math
import matplotlib as plt
import argparse

# define functions

def calculateContacts():
# needs contacts per frame, contacts per residue

def calculateTimeConstants():
# needs distance in z, time

def calculateRDF():

def calculateLipidDistances():

def calculatePowerFit():

def calculateLipidComposition():

def calculateClusterBindCorr():

def plotHistogram():

def plotLine():

# load and check files

traj = md.load('./*.xtc', top='./*.pdb') 

print "Shape of Cartesian coordinate array:"
print traj.xyz.shape
print "Time of final frame (ns):"
print traj.time[-1]

# save coordinates of protein + lipids + ions only:

atoms_to_keep = [a.index for a in traj.topology.atoms if a.name != 'W']
traj.restrict_atoms(atoms_to_keep)
traj.save('trajectory_for_analysis.h5')

# run analysis

## 1. Protein binding - time constants, contacts, effect on clustering

### time constants

### contacts

## 2. Lipid clustering  

### RDF

### clustering over time

## 3. Membrane curvature - fitting of power laws

### fit power law

## 4. correlations

### composition in defined regions of curvature

### clustering & binding

# produce plots
