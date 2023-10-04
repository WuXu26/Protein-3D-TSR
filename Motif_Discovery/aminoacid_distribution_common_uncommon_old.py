import glob, os, csv, ntpath,socket,argparse, time, re, gc
import pandas as pd, numpy as np
from collections import Counter
from joblib import Parallel, delayed, cpu_count
from os.path import expanduser
import itertools
import matplotlib as plt
import resource

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
parser.add_argument('--keys', '-keys', metavar='keys', \
	default='6362130',\
	help='Key that is to be searched.')

def search_by_key(fileName, file, keys, iskeys_file, outFolder):
	search_records = []
	not_search_records = []
	if iskeys_file == 1:
		line_count = len(file)
		df_com = file[file['key'].isin(keys)]
		un_df = file[~(file['key'].isin(keys))]
		s = time.time()
		df_com.to_csv(os.path.join(outFolder, fileName +'_common_keys.csv'))
		print('time writing:', (time.time() - s)/60)
		un_df.to_csv(os.path.join(outFolder, fileName +'_uncommon_keys.csv'))	
	else:
		common_out_file = open(os.path.join(outFolder, \
				fileName +'_common_keys.csv'),'w')
		uncommon_out_file = open(os.path.join(outFolder, \
					fileName +'_common_keys.csv'),'w')
		line_count = 0
		for line in file:
			line_count += 1
			line_splitted = line.split('\t')
			if line_splitted[0].strip() in keys:
				common_out_file.writelines(line_splitted)
				#search_records.append(line_splitted)
			else:
				uncommon_out_file.writelines(line_splitted)
				#not_search_records.append(line_splitted)
		common_out_file.close()
		uncommon_out_file.close()
		
		
	#return (fileName, line_count, len(set(df['key'].values)), len(df['key'].values))

def common_keys(files):
	print( 'Common Keys Calculation started.')
	common_keys = []
	start = time.time()
	for file in files:
		print( ntpath.basename(file)[:4])
		keys = []
		for line in open(file, 'r'):
			keys.append(line.split('\t')[0])
		if common_keys:
			common_keys = list(set(common_keys) & set(keys))
		else:
			common_keys = list(set(keys))
		
	print( 'Time taken for Common Keys Calculation: ', (time.time() - start)/60)
	return common_keys

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
def calculate_statistics(file, is_freq_groups, fileType):
		print('reading key file', file)
		use_col_names = ['key', 'aa0', 'aa1', 'aa2', 'theta', 'distance']
		if fileType == 'full':
			df_key = pd.read_table(file, delimiter = '\t', names = ['key', 'freq'], \
				dtype = {'key': np.dtype(int), 'freq':  np.dtype(int)})
			s = time.time()
			df_triplets = pd.read_csv(os.path.join(args.path, \
				args.sample_name, args.setting,\
				'{}.triplets_{}'.format(ntpath.basename(file)[:4], args.setting)) , \
				usecols = use_col_names, \
				dtype = {'key': np.dtype(int), 'aa0': str,  'aa1': str,  'aa2': str,\
				'theta': np.dtype(float), 'distance': np.dtype(float)}, \
				delimiter = '\t',names = column_names)
			print('reading triplet file(mins):', (time.time() - s)/60, \
				'mem:',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, )
			df = df_key.merge(df_triplets,on = 'key',how='left')
			del df_key, df_triplets
			if is_freq_groups == 0:
				df['freqGroups'] = df['freq'].apply(generateFrequencyGroups)
			else:
				df['freqGroups'] = df['distance'].apply(generateDistanceGroups)
			df.to_csv(os.path.join(args.path, args.sample_name, args.setting,\
				'{}.merged'.format(ntpath.basename(file)[:4])))
		else:
			df = pd.read_csv(file , delimiter = ',', header = 0)
		
		all_mean_theta = df['theta'].mean()
		all_mean_dist = df['distance'].mean()
		all_std_theta = df['theta'].std()
		all_std_dist = df['distance'].std()		
		
		#freqGroups_maxDist_avg = df
		#freqGroups_maxDist_avg = freqGroups_maxDist_avg[['distance','freqGroups']].groupby(['freqGroups']).agg('mean')
		#print(freqGroups_maxDist_avg)
		

		#freqGroups_theta_avg = df
		#freqGroups_theta_avg = freqGroups_theta_avg[['theta','freqGroups']].groupby(['freqGroups']).agg('mean')
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
		#aa_distribution_no_group = aminoacid_freq_no_grouping[['aacd', 'percent']].set_index('aacd')


		# aacd0 = df[['freq','aa0','freqGroups']]
		# aacd0.columns= ['freq','aacd','freqGroups']
		# aacd1 = df[['freq','aa1','freqGroups']]
		# aacd1.columns= ['freq','aacd','freqGroups']
		# aacd2 = df[['freq','aa2','freqGroups']]
		# aacd2.columns= ['freq','aacd','freqGroups']
		# groubyColumns = ['aacd','freqGroups']

		# aminoacid_freq = aacd0.append(aacd1).append(aacd2)
		# aminoacid_freq = aminoacid_freq.groupby(groubyColumns)['freq'].agg('sum').reset_index(name='sum')

		# file1_dict = aminoacid_freq[['sum','freqGroups']].groupby(['freqGroups']).agg('sum').to_dict()
		# aminoacid_freq['totals'] = aminoacid_freq['freqGroups'].apply(lambda x: file1_dict['sum'][x])
		# aminoacid_freq['percent'] = 100*aminoacid_freq['sum']/aminoacid_freq['totals']

		# return freqGroups_maxDist_avg.reset_index(), \
		# 	freqGroups_theta_avg.reset_index(), \
		# 	aminoacid_freq, aminoacid_freq_no_grouping, \
		# 	all_mean_theta, all_mean_dist, all_std_theta, all_std_dist

		return aminoacid_freq_no_grouping, \
			all_mean_theta, all_mean_dist, all_std_theta, all_std_dist

def get_aa_percents(keys_files, fileType):
	all_proteins_results = []
	for file in keys_files:
		print(ntpath.basename(file)[:4])
		if int(args.is_freq_groups) == 0:
			writer = pd.ExcelWriter(os.path.join(outFolder,\
			'{}_{}_freq_grouping.xlsx'.format(ntpath.basename(file)[:4], fileType)), \
			engine='xlsxwriter')
		else:
			writer = pd.ExcelWriter(os.path.join(outFolder,\
			'{}_{}_distance_grouping.xlsx'.format(ntpath.basename(file)[:4], fileType)), \
			engine='xlsxwriter')

		# freqGroups_maxDist_avg, \
		# freqGroups_theta_avg, \
		# aminoacid_freq, \
		# aminoacid_freq_no_grouping, \
		# all_mean_theta, all_mean_dist, all_std_theta, all_std_dist= calculate_statistics(file, int(args.is_freq_groups), fileType)


		aminoacid_freq_no_grouping, \
		all_mean_theta, all_mean_dist, all_std_theta, all_std_dist= calculate_statistics(file, int(args.is_freq_groups), fileType)

		# all_proteins_results.append(dict(zip(aminoacid_freq_no_grouping.aacd, aminoacid_freq_no_grouping.percent)).\
		# 	update({'all_mean_theta': all_mean_theta, 'all_mean_dist': all_mean_dist, \
		# 		'all_std_theta':all_std_theta, 'all_std_dist': all_std_dist, \
		# 		'protein': ntpath.basename(file)[:4]}))
		

		d1 = dict(zip(aminoacid_freq_no_grouping.aacd, aminoacid_freq_no_grouping.percent))
		
		all_proteins_results.append(dict(d1, **{'all_mean_theta': all_mean_theta, 'all_mean_dist': all_mean_dist, \
				'all_std_theta':all_std_theta, 'all_std_dist': all_std_dist, \
				'protein': ntpath.basename(file)[:4]}))

		# freqGroups_maxDist_avg.to_excel(writer,sheet_name='maxDist_distribution')
		# freqGroups_theta_avg.to_excel(writer,sheet_name='theta_distribution')
		# aminoacid_freq.to_excel(writer,sheet_name='aa_distribution')
		aminoacid_freq_no_grouping.to_excel(writer,sheet_name='all_aa_distribution')

		writer.close()
		# del freqGroups_maxDist_avg, \
		# freqGroups_theta_avg, \
		# aminoacid_freq, \
		# aminoacid_freq_no_grouping, \
		# all_mean_theta, all_mean_dist, all_std_theta, all_std_dist

		del aminoacid_freq_no_grouping, \
		all_mean_theta, all_mean_dist, all_std_theta, all_std_dist
	#print all_proteins_results
	return all_proteins_results

if __name__ == '__main__':
	gc.collect()
	start_time = time.time()
	args = parser.parse_args()

	iskeys_file = 1
	triplet_files = glob.glob(os.path.join(args.path,\
					 args.sample_name, args.setting, '*.triplets*'))
	column_names = ['key', 'aa0', 'pos0', 'aa1', 'pos1', 'aa2', 'pos2', \
	'classT', 'theta', 'classL', 'distance', 'x0', 'y0', 'z0', 'x1', 'y1',\
	'z1', 'x2', 'y2', 'z2']
	use_col_names = ['key', 'freq', 'aa0', 'aa1', 'aa2', 'theta', 'distance', 'freqGroups']

	if int(args.is_freq_groups) == 0:
		print("FREQUENCY GROUPING...")
		outFolder = os.path.join(args.path, args.sample_name, args.setting, \
						'Frequency_groups')
	else:
		print("DISTANCE GROUPING...")
		outFolder = os.path.join(args.path, args.sample_name, args.setting, \
						'Distance_groups')
	if not os.path.exists(outFolder):
			os.makedirs(outFolder)
	
	all_keys_files = glob.glob(os.path.join(args.path, args.sample_name, args.setting, '*.keys*'))

	#Getting AA Percent for full Triplet Files
	print("Getting AA Percent for full Triplet Files..")
	writer = pd.ExcelWriter(os.path.join(outFolder,\
			'aa_percent_distribution_all_proteins.xlsx'), engine='xlsxwriter')
	pd.DataFrame(get_aa_percents(all_keys_files, 'full')).to_excel(writer, sheet_name= 'all_triplets')
	print("Done.")

	#----------------------------------------------------------------------------------------------------
	#Getting AA Percent for Common Keys Files
	print("Getting AA Percent for Common Keys Files..")
	keys = common_keys(all_keys_files)
	print('Total common keys are: {}'.format(len(keys)))
	outFolder = os.path.join(args.path, args.sample_name, args.setting, \
				'aa_percent_files')
	if not os.path.exists(outFolder):
			os.makedirs(outFolder)
	pd.DataFrame(keys).to_csv(os.path.join(outFolder,'common_keys.csv'))

	print("Preparing common keys files for every protein..")
	merged_files = glob.glob(os.path.join(args.path, args.sample_name, args.setting, '*.merged*'))
	for file in merged_files:
		print(ntpath.basename(file)[:4])	
		#If reading with pandas
		if iskeys_file == 1:
			start = time.time()
			
			df = pd.read_table(file, delimiter = ',', header = 0,  \
				usecols = use_col_names, \
				dtype = {'key': np.dtype(int), 'freq': np.dtype(int), 'aa0': str,  \
				'aa1': str,  'aa2': str, 'freqGroups': str, 'theta': np.dtype(float), 'distance': np.dtype(float)})
			
			print('mem usage:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, \
				'len:', len(df), 'time reading:', (time.time() - start)/60)
			search_by_key(ntpath.basename(file)[:4], \
						df,\
						keys,1, outFolder)
		#line by line file read
		else:
			search_by_key(ntpath.basename(file)[:4],\
						open(file,'r'),\
						keys, 0, outFolder)
	keys_files = glob.glob(os.path.join(outFolder, '*_common_keys.csv'))
	pd.DataFrame(get_aa_percents(keys_files, 'common')).to_excel(writer, sheet_name= 'common')
	print("Done.")

	#----------------------------------------------------------------------------------------------------
	#Getting AA Percent for UnCommon Keys Files
	print("Getting AA Percent for UnCommon Keys Files..")
	keys_files = glob.glob(os.path.join(outFolder, '*_uncommon_keys.csv'))
	pd.DataFrame(get_aa_percents(keys_files, 'uncommon')).to_excel(writer, sheet_name= 'uncommon')
	print("Done.")

	
	end_time=time.time()
	total_time=((end_time)-(start_time))
	print("Task completed. Total Time taken(mins): {}".format(total_time/ 60))