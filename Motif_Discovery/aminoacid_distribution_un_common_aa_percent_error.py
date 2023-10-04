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
	"""Searches the file for the keys that are passed."""
	if os.path.exists(os.path.join(outFolder,'{}_common_keys.csv'.format(fileName))) and \
		os.path.exists(os.path.join(outFolder,'{}_uncommon_keys.csv'.format(fileName))):
		print('search_by_key() : Skipping as {}_common_key.csv and {}_uncommon_key.csv exists. '.\
			format(ntpath.basename(file)[:4]))
		return

	print('search_by_key() : {}'.format(fileName))
	search_records = []
	not_search_records = []
	if iskeys_file == 1:
		line_count = len(file)
		df_com = file[file['key'].isin(keys)]
		un_df = file[~(file['key'].isin(keys))]
		df_com.to_csv(os.path.join(outFolder, fileName +'_common_keys.csv'))
		un_df.to_csv(os.path.join(outFolder, fileName +'_uncommon_keys.csv'))	
	else:
		common_out_file = open(os.path.join(outFolder, \
				fileName +'_common_keys.csv'),'w')
		uncommon_out_file = open(os.path.join(outFolder, \
					fileName +'_uncommon_keys.csv'),'w')
		line_count = 0
		for line in file:
			line_count += 1
			line_splitted = line.split('\t')
			if line_splitted[0].strip() in keys:
				common_out_file.writelines(line_splitted)
			else:
				uncommon_out_file.writelines(line_splitted)
		common_out_file.close()
		uncommon_out_file.close()		

def common_keys(files):
	"""Calculates common keys among all key files in the speified folder."""
	print( 'common_keys() : Common Keys Calculation started.')
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

	print('common_keys: Total Common keys among files: {}'.format(len(common_keys)))	
	print( 'common_keys: Time taken for Common Keys Calculation: {}'.\
		format((time.time() - start)/60) )
	return common_keys

def thetaClass_( binBoundaries, value):
	"""Identifies and assigns bins for the boundaries."""
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
	"""Generates frequency bins with the boundaries passed as argument."""
	x = float(x)
	binBounds = list(map(float, args.freqBounds.split(',')))	
	return thetaClass_( binBounds, x)

def generateDistanceGroups(x):
	"""Generates distance bins with the boundaries passed as argument."""
	x = float(x)
	binBounds = list(map(float, args.distBounds.split(',')))
	return thetaClass_( binBounds, x)

def merge_key_triplet_file(key_file, triplet_file, file_name, is_freq_groups, delim, is_header):
	"""Merges keys file and triplet file by join on key number.
	Also add grouping details."""
	s = time.time()
	use_col_names = ['key', 'aa0', 'aa1', 'aa2', 'theta', 'distance']
	if not key_file:
		key_file = os.path.join(args.path, \
				args.sample_name, args.setting,\
				'{}.keys_{}'.format(file_name, args.setting))
	if not triplet_file:
		triplet_file = os.path.join(args.path, \
				args.sample_name, args.setting,\
				'{}.triplets_{}'.format(file_name, args.setting))

	df_key = pd.read_table(key_file, delimiter = '\t', names = ['key', 'freq'], \
				dtype = {'key': np.dtype(int), 'freq':  np.dtype(int)})
	
	if is_header == 0:
		df_triplets = pd.read_csv(triplet_file , \
					usecols = use_col_names, \
					dtype = {'key': np.dtype(int), 'aa0': str,  'aa1': str,  'aa2': str,\
					'theta': np.dtype(float), 'distance': np.dtype(float)}, \
					delimiter = delim, names = ['key', 'aa0', 'pos0', 'aa1', 'pos1', 'aa2', 'pos2', \
					'classT', 'theta', 'classL', 'distance', 'x0', 'y0', 'z0', 'x1', 'y1',\
					'z1', 'x2', 'y2', 'z2'])
	else:
		df_triplets = pd.read_csv(triplet_file , \
					usecols = use_col_names, \
					dtype = {'key': np.dtype(int), 'aa0': str,  'aa1': str,  'aa2': str,\
					'theta': np.dtype(float), 'distance': np.dtype(float)}, \
					delimiter = delim, header = 0)
	
	df = df_triplets.merge(df_key,on = 'key',how='left')
	del df_key, df_triplets
	if is_freq_groups == 0:
		df['freqGroups'] = df['freq'].apply(generateFrequencyGroups)
	else:
		df['freqGroups'] = df['distance'].apply(generateDistanceGroups)
	# print('merge_key_triplet_file() : {}- Time taken for merging key and triplet(mins): {}, mem: {}'.format(file_name, \
	# 	round((time.time() - s)/60, 2), \
	# 	resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, ))
	return df

def calculate_statistics(file, is_freq_groups, fileType):
		"""Calculates amino acid percentanges based on grouping."""
		#print('calculate_statistics() : Working on file', file)
		if fileType == 'full':
			if os.path.exists(os.path.join(args.path, args.sample_name, args.setting,\
				'{}.merged'.format(ntpath.basename(file)[:4]))):
				print('Merged file exists')
				df = pd. read_csv(os.path.join(args.path, args.sample_name, args.setting,\
				'{}.merged'.format(ntpath.basename(file)[:4])), delimiter = ',', header = 0)
			else:
				df = merge_key_triplet_file(file, None, ntpath.basename(file)[:4], is_freq_groups, '\t', 0)	
				df.to_csv(os.path.join(args.path, args.sample_name, args.setting,\
					'{}.merged'.format(ntpath.basename(file)[:4])))
		else:
			df = pd.read_csv(file , delimiter = ',', header = 0)
		
		all_mean_theta = df['theta'].mean()
		all_mean_dist = df['distance'].mean()
		all_std_theta = df['theta'].std()
		all_std_dist = df['distance'].std()		
		all_mean_freq = df['freq'].mean()
		all_std_freq = df['freq'].std()
				
		aacd0_no_grouping = df[['freq','aa0']]
		aacd0_no_grouping.columns= ['freq','aacd']
		aacd1_no_grouping = df[['freq','aa1']]
		aacd1_no_grouping.columns= ['freq','aacd']
		aacd2_no_grouping = df[['freq','aa2']]
		aacd2_no_grouping.columns= ['freq','aacd']
		aminoacid_freq_no_grouping = aacd0_no_grouping.append(aacd1_no_grouping).append(aacd2_no_grouping)
		aminoacid_freq_no_grouping = aminoacid_freq_no_grouping.groupby(['aacd'])['freq'].agg('sum').reset_index(name='sum')
		aminoacid_freq_no_grouping['percent'] = 100* aminoacid_freq_no_grouping['sum']/aminoacid_freq_no_grouping['sum'].sum()

		return aminoacid_freq_no_grouping, \
			all_mean_theta, all_mean_dist, all_std_theta, all_std_dist, all_mean_freq, all_std_freq

def get_aa_percents(keys_files, fileType):
	"""Get amino acid percentages for all files in the specified group."""
	all_proteins_results = []
	for file in keys_files:
		print('get_aa_percents() : Working on file {}'.format(ntpath.basename(file)[:4]))
		if int(args.is_freq_groups) == 0:
			writer = pd.ExcelWriter(os.path.join(outFolder,\
			'{}_{}_freq_grouping.xlsx'.format(ntpath.basename(file)[:4], fileType)), \
			engine='xlsxwriter')
		else:
			writer = pd.ExcelWriter(os.path.join(outFolder,\
			'{}_{}_distance_grouping.xlsx'.format(ntpath.basename(file)[:4], fileType)), \
			engine='xlsxwriter')

		aminoacid_freq_no_grouping, \
		all_mean_theta, all_mean_dist, all_std_theta, all_std_dist, \
		all_mean_freq, all_std_freq = calculate_statistics(file, int(args.is_freq_groups), fileType)
		
		d1 = dict(zip(aminoacid_freq_no_grouping.aacd, aminoacid_freq_no_grouping.percent))
		
		all_proteins_results.append(dict(d1, **{'all_mean_theta': all_mean_theta, 'all_mean_dist': all_mean_dist, \
				'all_std_theta':all_std_theta, 'all_std_dist': all_std_dist, 'all_mean_freq': all_mean_freq, 'all_std_freq':all_std_freq,\
				'protein': ntpath.basename(file)[:4]}))

		aminoacid_freq_no_grouping.to_excel(writer,sheet_name='all_aa_distribution')
		writer.close()

		del aminoacid_freq_no_grouping, \
		all_mean_theta, all_mean_dist, all_std_theta, all_std_dist
	return all_proteins_results

def generate_common_key_files(all_keys_files, outFolder):
	"""Generates common keys files under 'aa_percent_files' required for aminoacid percentage."""
	if os.path.exists(os.path.join(args.path, args.sample_name, args.setting, \
						'common_keys')):
		print('generate_common_key_files() : Common Key files already exist from a previous program. Adding frequency and grouping to the common key files.')
		common_triplet_files = glob.glob(os.path.join(os.path.join(args.path, args.sample_name, args.setting, \
						'common_keys'), '*_key_search_common_keys.csv'))
		if(len(common_triplet_files) != len(all_keys_files)):
			print("ERROR: generate_common_key_files() : All proteins did not have common keys triplet files generated. Repeat common key files calculation. Aborting..")
			return
		
		for file in common_triplet_files:
			if os.path.exists(os.path.join(outFolder,\
				'{}_common_keys.csv'.format(ntpath.basename(file)[:4]))):
				print('generate_common_key_files() : Skipping as {}_common_key.csv exists. '.format(ntpath.basename(file)[:4]))
				continue
			df = merge_key_triplet_file(None, file, ntpath.basename(file)[:4], int(args.is_freq_groups), ',', 1)
			df.to_csv(os.path.join(outFolder,\
				'{}_common_keys.csv'.format(ntpath.basename(file)[:4])))

	else:
		print('generate_common_key_files() : Common Key files not present. Generating common key files.')
		use_col_names = ['key', 'freq', 'aa0', 'aa1', 'aa2', 'theta', 'distance', 'freqGroups']
		keys = common_keys(all_keys_files)
		pd.DataFrame(keys).to_csv(os.path.join(outFolder,'common_keys.csv'))

		merged_files = glob.glob(os.path.join(args.path, args.sample_name, args.setting, '*.merged*'))
		for file in merged_files:
			#If reading with pandas
			if iskeys_file == 1:
				start = time.time()
				
				df = pd.read_table(file, delimiter = ',', header = 0,  \
					usecols = use_col_names, \
					dtype = {'key': np.dtype(int), 'freq': np.dtype(int), 'aa0': str,  \
					'aa1': str,  'aa2': str, 'freqGroups': str, 'theta': np.dtype(float), 'distance': np.dtype(float)})		
				# print('mem usage:', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, \
				# 	'len:', len(df), 'time reading:', (time.time() - start)/60)
				search_by_key(ntpath.basename(file)[:4], \
							df,\
							keys,1, outFolder)
			#line by line file read
			else:
				search_by_key(ntpath.basename(file)[:4],\
							open(file,'r'),\
							keys, 0, outFolder)

def generate_uncommon_key_files(all_keys_files, outFolder, is_freq_groups):
	"""Generates common keys files under 'aa_percent_files' required for aminoacid percentage."""
	setting_folder = os.path.join(args.path,\
					 args.sample_name, args.setting)
	un_common_key_files = glob.glob(os.path.join(outFolder, '*_uncommon_keys.csv'))
	if not un_common_key_files:
		print('generate_uncommon_key_files() : Generating Uncommon key files')
		for file in all_keys_files:
			print file
			file_name = ntpath.basename(file)[:4]
			c_file = os.path.join(outFolder, '{}_common_keys.csv'.format(file_name))
			m_file = os.path.join(setting_folder, '{}.merged'.format(file_name))
			if not c_file:
				print('generate_uncommon_key_files() : {} doesnot have common key file generated. Aborting..'.format(file_name))
				return
			c_keys = list(pd.read_table(c_file, delimiter = ',', header = 0,  \
					dtype = {'key': np.dtype(int), 'freq': np.dtype(int), 'aa0': str,  \
					'aa1': str,  'aa2': str, 'freqGroups': str, 'theta': np.dtype(float), \
					'distance': np.dtype(float)}).key.unique())
			#c_keys = list(common_key_file.key.unique())
			if os.path.exists(os.path.join(setting_folder, '{}.merged'.format(file_name))):
				print('generate_uncommon_key_files() : {}  Subtract from .merged file'.format(file_name))
				m_file = pd.read_table(m_file, delimiter = ',', header = 0,  \
					dtype = {'key': np.dtype(int), 'freq': np.dtype(int), 'aa0': str,  \
					'aa1': str,  'aa2': str, 'freqGroups': str, 'theta': np.dtype(float), \
					'distance': np.dtype(float)})
			else:
				print('generate_uncommon_key_files() : {} Subtract from .triplets file'.format(file_name))
				m_file = merge_key_triplet_file(file, None, file_name, is_freq_groups, '\t', 0)
			un_df = m_file[~(m_file['key'].isin(c_keys))]
			un_df.to_csv(os.path.join(outFolder, file_name +'_uncommon_keys.csv'))

	else:
		print('generate_uncommon_key_files() : Uncommon key files already exist.')
		return

if __name__ == '__main__':
	gc.collect()
	start_time = time.time()
	args = parser.parse_args()

	iskeys_file = 1
	triplet_files = glob.glob(os.path.join(args.path,\
					 args.sample_name, args.setting, '*.triplets*'))	

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

	#----------------------------------------------------------------------------------------------------
	#Getting AA Percent for full Triplet Files
	print('main(): #----------------------------------------------------------------------------------------------------')
	print("main() : Getting AA Percent for full Triplet Files..")
	writer = pd.ExcelWriter(os.path.join(outFolder,\
			'aa_percent_distribution_all_proteins.xlsx'), engine='xlsxwriter')
	pd.DataFrame(get_aa_percents(all_keys_files, 'full')).to_excel(writer, sheet_name= 'all_triplets')#Uncomment this:Sarika
	print("main() : Done. Getting AA Percent for full Triplet Files.")
	print('main(): #----------------------------------------------------------------------------------------------------')

	#----------------------------------------------------------------------------------------------------
	outFolder = os.path.join(args.path, args.sample_name, args.setting, \
					'aa_percent_files')
	if not os.path.exists(outFolder):
				os.makedirs(outFolder)
	#----------------------------------------------------------------------------------------------------			
	#Getting AA Percent for Common Keys Files
	print('main(): #----------------------------------------------------------------------------------------------------')
	print("main() : Getting AA Percent for Common Keys Files..")	
	generate_common_key_files(all_keys_files, outFolder)

	keys_files = glob.glob(os.path.join(outFolder, '*_common_keys.csv'))
	pd.DataFrame(get_aa_percents(keys_files, 'common')).to_excel(writer, sheet_name= 'common')
	print("main() : Done. Getting AA Percent for Common Keys Files.")
	print('main(): #----------------------------------------------------------------------------------------------------')

	#----------------------------------------------------------------------------------------------------
	#Getting AA Percent for UnCommon Keys Files
	print('main(): #----------------------------------------------------------------------------------------------------')
	print("main() : Getting AA Percent for UnCommon Keys Files..")
	generate_uncommon_key_files(all_keys_files, outFolder, args.is_freq_groups)
	keys_files = glob.glob(os.path.join(outFolder, '*_uncommon_keys.csv'))
	pd.DataFrame(get_aa_percents(keys_files, 'uncommon')).to_excel(writer, sheet_name= 'uncommon')
	print("main() : Done. Getting AA Percent for UnCommon Keys Files.")
	print('main(): #----------------------------------------------------------------------------------------------------')

	
	end_time=time.time()
	total_time=((end_time)-(start_time))
	print("main() : Task completed. Total Time taken(mins): {}".format(total_time/ 60))