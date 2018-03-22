# import generic python modules

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

import mdtraj as mdt
import numpy as np
import subprocess
import math
import matplotlib as plt
import argparse

# options

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

# store inputs

# parsing
args = parser.parse_args()

args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]

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
