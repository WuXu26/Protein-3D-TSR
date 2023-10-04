import glob, os, csv, ntpath,socket,argparse, time, re
import pandas as pd, numpy as np
from collections import Counter
from joblib import Parallel, delayed, cpu_count
from os.path import expanduser
import itertools
import matplotlib as plt

__author__ = "Venkata Sarika Kondra"
__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"

parser = argparse.ArgumentParser(description='Search Keys.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', \
	default='t1', \
	help='Name of the sample on which this script should be run.')
parser.add_argument('--path', '-path', metavar='path', \
	default=os.path.join(expanduser('~'),'Research', 'Protien_Database', \
		'extracted_new_samples', 'testing'), \
	help='Directory of input sample and other files.')
parser.add_argument('--keys', '-keys', metavar='keys', \
	default='6362130',\
	help='Key that is to be searched.')
parser.add_argument('--setting', '-setting', metavar='setting', \
	default='theta29_dist35', \
	help='Setting with theta and distance bin values.')
parser.add_argument('--is_freq_groups', '-is_freq_groups', metavar='is_freq_groups', \
	default=0, \
	help='0 = Distribution w.r.t frequency of key, 1 = Distribution w.r.t  max distance.')
parser.add_argument('--freqBounds', '-freqBounds', metavar='freqBounds', \
    default = '5,200', \
    help='Bin Boundaries for frequency.')
parser.add_argument('--distBounds', '-distBounds', metavar='distBounds', \
    default = '7,9,11,15,20', \
    help='Bin Boundaries for maxDist.')

def thetaClass_( binBoundaries, value):
	classL = -1
	if value >= binBoundaries[-1]: 
		return '>{}'.format(int(binBoundaries[-1]))
	if value < binBoundaries[0]:#Bins are seperately handled for theta and maxdist.
		return '<{}'.format(int(binBoundaries[0]))
	for i in binBoundaries:
		if (value < i) :#If the value is less than the boundary it falls in previous bin.
			classL = '{}-{}'.format(int(binBoundaries[binBoundaries.index(i) - 1]), int(i))
			break	
	return classL
def generateFrequencyGroups(x):
	x = float(x)
	binBounds = list(map(float, args.freqBounds.split(',')))	
	return thetaClass_( binBounds, x)
	# group = ''
	# if x <5:
	# 	group = '<5'
	# elif 5 <= x < 200:
	# 	group = '6-200'
	# else:
	# 	group = '>200'
	#return group

def generateDistanceGroups(x):
	x = float(x)
	binBounds = list(map(float, args.distBounds.split(',')))
	return thetaClass_( binBounds, x)	
	# group = ''
	# if x <7:
	# 	group = 'short'
	# elif 7 <= x < 15:
	# 	group = 'medium'
	# else:
	# 	group = 'long'
	# return group
def calculate_statistics(file, is_freq_groups):
		df_key = pd.read_table(file, delimiter = '\t', names = ['key', 'freq'])
		
		df_triplets = pd.read_table(os.path.join(args.path, \
			args.sample_name, \
			args.setting,\
			'{}.triplets_{}'.format(ntpath.basename(file)[:4], args.setting)) , delimiter = '\t',names = column_names)
		df_triplets = df_triplets[['key','aa0', 'pos0', 'aa1', 'pos1', 'aa2', 'pos2','classT','theta','classL','distance']]
		df = df_key.merge(df_triplets,on = 'key',how='left')
		if is_freq_groups == 0:
			df['freqGroups'] = df['freq'].apply(generateFrequencyGroups)
		else:
			df['freqGroups'] = df['distance'].apply(generateDistanceGroups)
		df.to_csv(os.path.join(args.path, \
			args.sample_name, \
			args.setting,\
			'{}.merged'.format(ntpath.basename(file)[:4])))
		freqGroups_maxDist_avg = df
		freqGroups_maxDist_avg = freqGroups_maxDist_avg[['distance','freqGroups']].groupby(['freqGroups']).agg('mean')
		#print(freqGroups_maxDist_avg)
		

		freqGroups_theta_avg = df
		freqGroups_theta_avg = freqGroups_theta_avg[['theta','freqGroups']].groupby(['freqGroups']).agg('mean')
		#print(freqGroups_theta_avg)
		
		aacd0_no_grouping = df[['freq','aa0']]
		aacd0_no_grouping.columns= ['freq','aacd']
		aacd1_no_grouping = df[['freq','aa1']]
		aacd1_no_grouping.columns= ['freq','aacd']
		aacd2_no_grouping = df[['freq','aa2']]
		aacd2_no_grouping.columns= ['freq','aacd']
		aminoacid_freq_no_grouping = aacd0_no_grouping.append(aacd1_no_grouping).append(aacd2_no_grouping)
		aminoacid_freq_no_grouping = aminoacid_freq_no_grouping.groupby(['aacd'])['freq'].agg('sum').reset_index(name='sum')
		aminoacid_freq_no_grouping['percent'] = 100* aminoacid_freq_no_grouping['sum']/aminoacid_freq_no_grouping['sum'].sum()
		#print aminoacid_freq_no_grouping

		aacd0 = df[['freq','aa0','freqGroups']]
		aacd0.columns= ['freq','aacd','freqGroups']
		aacd1 = df[['freq','aa1','freqGroups']]
		aacd1.columns= ['freq','aacd','freqGroups']
		aacd2 = df[['freq','aa2','freqGroups']]
		aacd2.columns= ['freq','aacd','freqGroups']
		groubyColumns = ['aacd','freqGroups']

		aminoacid_freq = aacd0.append(aacd1).append(aacd2)
		aminoacid_freq = aminoacid_freq.groupby(groubyColumns)['freq'].agg('sum').reset_index(name='sum')

		file1_dict = aminoacid_freq[['sum','freqGroups']].groupby(['freqGroups']).agg('sum').to_dict()
		aminoacid_freq['totals'] = aminoacid_freq['freqGroups'].apply(lambda x: file1_dict['sum'][x])
		aminoacid_freq['percent'] = 100*aminoacid_freq['sum']/aminoacid_freq['totals']

		#print aminoacid_freq		
		return freqGroups_maxDist_avg.reset_index(), \
			freqGroups_theta_avg.reset_index(), \
			aminoacid_freq, aminoacid_freq_no_grouping

if __name__ == '__main__':
	start = time.time()
	args = parser.parse_args()

	#binBounds = list(map(float, args.binBounds.split(',')))

	#print thetaClass_(binBounds, 208), thetaClass_(binBounds, 2) , thetaClass_(binBounds, 6) , thetaClass_(binBounds, 176)
	triplet_files = glob.glob(os.path.join(args.path,\
					 args.sample_name, args.setting, '*.triplets*'))
	column_names = ['key', 'aa0', 'pos0', 'aa1', 'pos1', 'aa2', 'pos2', \
	'classT', 'theta', 'classL', 'distance', 'x0', 'y0', 'z0', 'x1', 'y1',\
	'z1', 'x2', 'y2', 'z2']

	if int(args.is_freq_groups) == 0:
		outFolder = os.path.join(args.path, args.sample_name, args.setting, \
						'Frequency_groups')
	else:
		outFolder = os.path.join(args.path, args.sample_name, args.setting, \
						'Distance_groups')
	if not os.path.exists(outFolder):
			os.makedirs(outFolder)
	
	keys_files = glob.glob(os.path.join(args.path, args.sample_name, args.setting, '*.keys*'))

	for file in keys_files:
		print(ntpath.basename(file)[:4])
		if int(args.is_freq_groups) == 0:
			writer = pd.ExcelWriter(os.path.join(outFolder,\
			'{}_freq_grouping.xlsx'.format(ntpath.basename(file)[:4])), \
			engine='xlsxwriter')
		else:
			writer = pd.ExcelWriter(os.path.join(outFolder,\
			'{}_distance_grouping.xlsx'.format(ntpath.basename(file)[:4])), \
			engine='xlsxwriter')

		freqGroups_maxDist_avg, freqGroups_theta_avg, aminoacid_freq, aminoacid_freq_no_grouping = calculate_statistics(file, int(args.is_freq_groups))

		freqGroups_maxDist_avg.to_excel(writer,sheet_name='maxDist_distribution')
		freqGroups_theta_avg.to_excel(writer,sheet_name='theta_distribution')
		aminoacid_freq.to_excel(writer,sheet_name='aa_distribution')
		aminoacid_freq_no_grouping.to_excel(writer,sheet_name='all_aa_distribution')

		writer.close()
	end_time=time.time()
	total_time=((end_time)-(start_time))
	print("Task completed. Total Time taken(secs): {}".format(total_time))