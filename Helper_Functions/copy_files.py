#!/usr/bin/env python
#copy_files.py
""" 
	Copy protein files in Details file to required_files folder in same path.

"""

import shutil
import argparse
import os
import pandas as pd
from os.path import expanduser

parser = argparse.ArgumentParser(description='Copy files to folder.')
parser.add_argument('--path', '-path', metavar='path', \
	default=os.path.join(expanduser('~'),'Research', 'Protien_Database', \
		'extracted_new_samples', 'testing'), \
	help='Directory of input sample and other files.')

args = parser.parse_args()

outFolder = os.path.join(args.path, 'required_files')
if not os.path.exists(outFolder):
	os.makedirs(outFolder)

df = pd.read_csv(os.path.join(args.path, 'sample_details_short_list.csv'), header = 0)
for protein in df['protein'].values:
	if os.path.exists(os.path.join(args.path,"{}.pdb".format(protein))):
		shutil.copy(os.path.join(args.path,"{}.pdb".format(protein)), outFolder)

print("Copying files done to location {}".format(outFolder))