#!/usr/bin/env python

import itertools
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='convert lammps trajectory to xyz format')
parser.add_argument('-i','--input',help ="input lammps file to be read",default = 'coords.lammpstrj')
args = parser.parse_args()

filename = args.input

f_frames = open('frames.dat','w')
f_coords = open('coords.dat','w')

nlines = 9
frames = 1
write = False
with open(filename) as f:
    for line in f:
        if line.startswith('ITEM: TIMESTEP'):
            write = False
            f_frames.write(str(frames)+'\n')
            frames = frames+1

        if(write):
            line=line.strip()
            columns=line.split()
            txt=columns[1]+' '+columns[2]+' '+columns[3]+' '+columns[4]+'\n'
            f_coords.write(txt)

        if line.startswith('ITEM: ATOMS'):
            write = True

f_frames.close()
f_coords.close()
#natoms = int(natoms)
#print 'Number of atoms ',natoms
#print 'box size ', high-low
