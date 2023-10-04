import numpy as np
import os
import glob
import argparse
import ntpath


__author__ = "Venkata Sarika Kondra"

__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"


parser = argparse.ArgumentParser(description='Calculate shortest distance between keys.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', default='t', \
    help='Name of the sample on which this script should be run.')
parser.add_argument('--path', '-path', metavar='path', \
    default='/media/c00219805/Elements1/Research/Hierarchical_Classification/Database/', \
    help='Directory of input sample and other files.')
parser.add_argument('--bin_folder', '-bin_folder', metavar='bin_folder', \
    default='theta29_dist35', \
    help='Directory of input sample and other files.')

parser.add_argument('--threshold_dist', '-threshold_dist', metavar='threshold_dist', \
    default='10', \
    help='Threshold distance.')

parser.add_argument('--k1_line', '-k1_line', metavar='k1_line', \
    default='4590902	LYS	  10 	GLU	  12 	PHE	  11 	29	89.2906413036	2	5.58234950536	39.634	30.133	20.962	38.905	34.036	24.886	36.868	32.588	22.026', \
    help='full line f key under search.')

def calculate_distance(v_a, v_b):
    p1 = np.array(v_a)
    p2 = np.array(v_b)

    squared_dist = np.sum((p1-p2)**2, axis=0)
    dist = np.sqrt(squared_dist)
    return dist

def process_lines(line):
	return

if __name__ == '__main__':

	args = parser.parse_args()

	triplet_files = glob.glob(os.path.join(args.path, args.sample_name, 
		args.bin_folder, 
		'*.triplets_{}'.format(args.bin_folder)))
	print(triplet_files)

	values_k1 = args.k1_line.replace("\n","").split('\t')
	k1 = [(float(values_k1[11]), float(values_k1[12]), float(values_k1[13])), 
			(float(values_k1[14]), float(values_k1[15]), float(values_k1[16])), 
			(float(values_k1[17]), float(values_k1[18]), float(values_k1[19]))]


	for file in triplet_files:
		out_file = open(os.path.join(args.path, args.sample_name, args.bin_folder,'{}.{}_{}'.format(
			ntpath.basename(file).split('.')[0], args.threshold_dist, values_k1[0])), "w")
		with open(file) as f:
			for line in f.readlines():
				values = line.split('\t')
				k2 = [(float(values[11]), float(values[12]), float(values[13])), 
							(float(values[14]), float(values[15]), float(values[16])), 
							(float(values[17]), float(values[18]), float(values[19]))]
				short_dist = calculate_distance(k1[0],k2[0])
				for cord1 in k1:
					for cord2 in k2:
					    #print(calculate_distance(cord1,cord2))
					    short_dist = min(short_dist,calculate_distance(cord1,cord2))
				if short_dist > 0 and short_dist <= float(args.threshold_dist) :
					    	print(short_dist)
					        out_file.writelines (line.replace("\n", "\t") + "{}".format(short_dist) + "\n")
		out_file.close()
