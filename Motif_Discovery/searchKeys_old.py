#!/usr/bin/env python
#searchKeys.py
"""
	Searches for specific keys or amino acids or positions and retrieves 
	their details from triplets files
"""

import glob, os, csv, ntpath,socket,argparse, time, re
import pandas as pd, numpy as np
from collections import Counter
from joblib import Parallel, delayed, cpu_count
from os.path import expanduser
import itertools


__author__ = "Venkata Sarika Kondra"
__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"

parser = argparse.ArgumentParser(description='Search Keys.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', \
	default='sample_p_loop_10', \
	help='Name of the sample on which this script should be run.')
parser.add_argument('--path', '-path', metavar='path', \
	default=os.path.join(expanduser('~'),'Research', 'Protien_Database', \
		'extracted_new_samples', 'testing'), \
	help='Directory of input sample and other files.')
parser.add_argument('--keys', '-keys', metavar='keys', \
	default='6362130,6362129,6362128,6362131,6362132',\
	help='Key that is to be searched.')
parser.add_argument('--aacds', '-aacds', metavar='aacds', \
	default='gly,thr,lys',\
	help='Amino acids that are to be searched.')
parser.add_argument('--setting', '-setting', metavar='setting', \
	default='theta29_dist35', \
	help='Setting with theta and distance bin values.')
parser.add_argument('--exclude', '-exclude', metavar='exclude', \
	default=True, \
	help='Exclude certain keys from when extracting triplets details.')
parser.add_argument('--search_mode', '-search_mode', metavar='search_mode', \
	default=2, \
	help='0 if Key search, 1 if amino acid search and 2 if sequence identification.')

def search_by_key(file,keys):
	print('inside search key..filtering started..')
	search_records = []
	line_count = len(file)
	df = file[file['key'].isin(keys)]
	# for line in file:
	# 	line_count += 1
	# 	line_splitted=line.split('\t')
	# 	if line_splitted[0].strip() in keys:
	# 		search_records.append(line_splitted)				
	#return pd.DataFrame(search_records,columns = column_names)
	return df, line_count, len(set(df['key'].values))

def search_aacd(file,aacds, req_pos):
	search_records = []	
	pattern_line = []
	for line in file:
		line_aacds = {}
		line_splitted=line.split('\t')
		line_aacd_arr = [line_splitted[1].strip().upper(), \
			line_splitted[3].strip().upper(), line_splitted[5].strip().upper()]	
		if sorted(line_aacd_arr) == sorted(aacds):
			search_records.append(line_splitted)
			if req_pos:
				line_pos_arr = [line_splitted[2].strip().upper(), \
					line_splitted[4].strip().upper(), \
					line_splitted[6].strip().upper()]
				if all(elem in req_pos for elem in line_pos_arr):
					pattern_line.append(line)
			
	return pd.DataFrame(search_records,columns = column_names), pattern_line

def identify_pattern_by_regex(files,column_names, patterns):
	matched = []
	all_pattern = generate_sequences_from_triplets_by_FR(file)
	for pattern in patterns:
		print(re.findall(pattern, all_pattern))
		matched.append(pattern)
	return ntpath.basename(file)[:4], matched

def identify_pattern_by_normal(file, pattern,count_of_x):
	matched = []
	# If file size is greater than 5 GB pandas is throwing Memory error
	# But if you always read line-by-line it is very slow. 
	# Hence using two way reading
	if float(os.path.getsize(file)/ 2**30) > 4 :
		print('File Size too big')
		return ntpath.basename(file)[:4], ['File Size too big']
		#all_pattern = generate_sequences_from_triplets_by_FR(file)	
	all_pattern = generate_sequences_from_triplets_by_DF(file)
	i = 0
	while i < len(all_pattern):
		if (all_pattern[i].split('_')[0] == pattern[0]) and ( i + count_of_x + len(pattern) < len(all_pattern)):
			#print 'checking this: {}'.format(all_pattern[i:i + count_of_x + len(pattern)])
			if (all_pattern[i + count_of_x +1].split('_')[0] == pattern[1]) \
				and (all_pattern[i+count_of_x +2].split('_')[0] == pattern[2]) \
				and (all_pattern[i+count_of_x +3].split('_')[0] == pattern[3]):
				print( 'found', all_pattern[i:i + count_of_x + len(pattern)])
				matched.append(all_pattern[i:i + count_of_x + len(pattern)])
		i += 1
	return (ntpath.basename(file)[:4], matched)

def generate_sequences_from_triplets_by_DF(file):
	print(ntpath.basename(file)[:4])
	all_pattern = []
	df = pd.read_table(file, delimiter = '\t', names = column_names)
	a_0 = pd.Series(df.aa0.values,index=df.pos0).to_dict()
	a_1 = pd.Series(df.aa1.values,index=df.pos1).to_dict()
	a_2 = pd.Series(df.aa2.values,index=df.pos2).to_dict()
	a_0.update(a_1) 
	a_0.update(a_2) 
	pos_aa_dict = sorted(a_0.items())
	for key, value in pos_aa_dict:
		all_pattern .append(str(value) + '_' + str(key))
	return all_pattern
def generate_sequences_from_triplets_by_FR(file):
	print( ntpath.basename(file)[:4])
	all_pattern = []
	pos_aa = {}
	for line in open(file,'r'):
		line_splitted=line.split('\t')
		if line_splitted[2].strip() not in pos_aa.keys():
			pos_aa[line_splitted[2].strip()] = line_splitted[1].strip().upper()
		if line_splitted[4].strip() not in pos_aa.keys():
			pos_aa[line_splitted[4].strip()] = line_splitted[3].strip().upper()
		if line_splitted[6].strip() not in pos_aa.keys():
			pos_aa[line_splitted[6].strip()] = line_splitted[5].strip().upper()
	for key, value in sorted(pos_aa.items()):
		all_pattern .append(str(value) + '_' + str(key))
	return all_pattern


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

def calculate_low_less15_freq_from_common_keys():
	writer = pd.ExcelWriter(args.path + args.sample_name + args.setting + \
		'common_keys/'+'//all_common_key_distribution.xlsx', engine='xlsxwriter')
	writer_low_freq = pd.ExcelWriter(args.path + args.sample_name + args.setting + \
		'common_keys/'+'//only_low_key_distribution.xlsx', engine='xlsxwriter')
	writer_dist_less15 = pd.ExcelWriter(args.path + args.sample_name + args.setting + \
		'common_keys/'+'//dist_less15_triplets.xlsx', engine='xlsxwriter')
	summary = []
	commom_key_freqs = Counter()
	common_keys_files = glob.glob(args.path + args.sample_name + args.setting + \
		'common_keys/' +'*_key_search_common_keys.csv')
	for file in common_keys_files:
		print( ntpath.basename(file)[:4])
		df = pd.read_csv(file, delimiter = ',')
		x = Counter(df['key'].values)
		ckeys_freqs_df = pd.DataFrame(sorted(x.items(), key=lambda pair: pair[1], \
			 reverse=True), columns = ['key', 'freq'])
		ckeys_freqs_df.to_excel(writer,sheet_name=ntpath.basename(file)[:4])
		#Low Frequency Common Keys
		low_freqs = ckeys_freqs_df[ckeys_freqs_df['freq'] < 3]
		low_freqs.to_excel(writer_low_freq,sheet_name=ntpath.basename(file)[:4])
		df_low_freqs_triplets = df[ df['key'].isin(set(low_freqs['key'].values))]
		#Distance less than 15
		df_distance_less15 = df_low_freqs_triplets[df_low_freqs_triplets['distance'] <= 15]
		df_distance_less15.to_excel(writer_dist_less15,sheet_name=ntpath.basename(file)[:4])
		summary.append((ntpath.basename(file)[:4], len(set(low_freqs['key'].values)), \
			len(set(df_distance_less15['key'].values))))
		print (ntpath.basename(file)[:4], len(set(low_freqs['key'].values)), \
			len(set(df_distance_less15['key'].values)))
		#Overall Common Keys
		if commom_key_freqs:
			commom_key_freqs = commom_key_freqs + x
		else:
			commom_key_freqs = x

	pd.DataFrame(sorted(commom_key_freqs.items(), key=lambda pair: pair[1], reverse=True)). \
		to_excel(writer,sheet_name='Summary')

	df_merge = pd.read_csv(args.path + args.sample_name + args.setting +'common_keys/'+ \
		'summary_common_keys-7-31-2018.csv').merge(pd.DataFrame(summary, \
			columns = ['fileName', '# of low freq common keys', '# of low freq common keys with maxdist <= 15']),on = 'fileName')
	df_merge.to_csv(args.path + args.sample_name + args.setting + 'common_keys/' + 'summary_common_keys-8-02-2018.csv')

def process_files(file, req_aacds):
	print(ntpath.basename(file)[:4])
	pattern_theta_1 = ''
	pattern_dist_1 = ''
	pattern_theta_2 = ''
	pattern_dist_2 = ''
	
	samples_file = pd.read_csv(os.path.join(args.path, \
		args.sample_name, 'sample_details.csv'))
	req_pos = re.findall(r'\d+', samples_file.set_index('protein').\
		loc[ntpath.basename(file)[:4],  'pattern'])
	records, patterns = search_aacd(open(file,'r'),req_aacds, req_pos)
	if patterns:
		# pattern_theta_1 = [pattern.split('\t')[8] for pattern in patterns]
		# pattern_dist_1 = [pattern.split('\t')[10] for pattern in patterns]
		pattern_theta_1 = patterns[0].split('\t')[8]
		pattern_dist_1 = patterns[0].split('\t')[10]
		if len(patterns) >1:
			pattern_theta_2 = patterns[1].split('\t')[8]	
			pattern_dist_2 = patterns[1].split('\t')[10]
	records.to_csv(os.path.join(args.path, args.sample_name, \
			args.setting, 'common_keys', \
	 		ntpath.basename(file)[:4] + '_'+ "_".join(req_aacds) + '.csv'))
	print (ntpath.basename(file)[:4], pd.to_numeric(records['theta']).mean(),\
	 pd.to_numeric(records['distance']).mean(), pattern_theta_1, pattern_theta_2,\
	 pattern_dist_1, pattern_dist_2, patterns)
	return (ntpath.basename(file)[:4], pd.to_numeric(records['theta']).mean(),\
	 pd.to_numeric(records['distance']).mean(), pattern_theta_1, pattern_theta_2,\
	 pattern_dist_1, pattern_dist_2, patterns)

if __name__ == '__main__':
	start = time.time()
	args = parser.parse_args()
	files = glob.glob(os.path.join(args.path,\
					 args.sample_name, args.setting, '*.triplets*'))
	column_names = ['key', 'aa0', 'pos0', 'aa1', 'pos1', 'aa2', 'pos2', \
	'classT', 'theta', 'classL', 'distance', 'x0', 'y0', 'z0', 'x1', 'y1',\
	 'z1', 'x2', 'y2', 'z2']

	#files = glob.glob(args.path + args.sample_name + args.setting + '//*.keys*')
	#column_names = ['key','freq']

	#results_key_search = pd.DataFrame(columns = column_names)
	#summary = []	

	#Common Keys
	# keys = common_keys(files)

	# #Keys Search   
	if args.search_mode == 0:
		keys = [ str(key) for key in args.keys.split(',')]		    
		for file in files:
			print( ntpath.basename(file)[:4])
			start = time.time()
			print('inside common keys. Started reading file..')
			df_all = pd.read_table(file, delimiter = '\t', names = column_names)
			print('done reading file', (time.time() - start)/60)
			#records, line_count, matched_keys_count = search_key(open(file,'r'),keys)
			records, line_count, matched_keys_count = search_key(df_all,keys)
				
			records['fileName'] = ntpath.basename(file)[:4]
			summary.append([ntpath.basename(file)[:4], matched_keys_count, line_count])
			#results_key_search = results_key_search.append(records)	
			records.to_csv(args.path + args.sample_name + args.setting +'common_keys/'+ \
				ntpath.basename(file)[:4]+'_key_search_common_keys.csv')
		pd.DataFrame(summary, columns = ['fileName','common_keys','all_keys']) \
				.to_csv(args.path + args.sample_name + args.setting +'summary_common_keys.csv')

	#calculate_low_less15_freq_from_common_keys()

	# #Amino Acids Search
	if args.search_mode == 1:
		all_list = ['GLY', 'LYS','SER', 'GLY']
		already_completed = []
		for req_aacds in itertools.combinations(all_list, 3):
			if sorted(req_aacds) not in already_completed:
				already_completed.append(sorted(req_aacds))
				print('Started AminoAcid serach for: {}'.format("_".join(req_aacds)))
				means_theta = []
				means_distance = []
				writer = pd.ExcelWriter(os.path.join(args.path, args.sample_name, \
					args.setting, 'means_' + "_".join(req_aacds) + '.xlsx'), \
					engine = 'xlsxwriter')
				
				means_theta_distance  = Parallel(n_jobs=cpu_count() - 1, verbose=10, \
					backend="multiprocessing", batch_size="auto")(\
					delayed(process_files)(file, req_aacds) for file in files)
				pd.DataFrame(means_theta_distance, \
					columns = ['file', 'all_theta_per_protein', \
							'all_dist_per_protein', 'motif_theta_1','motif_theta_2', \
							'motif_dist_1', 'motif_dist_2','motif'])\
					.to_excel(writer,sheet_name='means_theta_distance')
				writer.close()

	#Sequence Identification
	if args.search_mode == 2:
		writer_pattern = pd.ExcelWriter(os.path.join(args.path, args.sample_name, \
					args.setting, 'patterns{}'.format('.xlsx')), \
					engine = 'xlsxwriter')
		patterns = [['GLY','GLY','LYS','SER'], ['GLY','GLY','LYS','THR']]
		for pattern in patterns:
			matched_list = []
			# sequences  = Parallel(n_jobs=cpu_count() - 1, verbose=10, \
			# 			backend="multiprocessing", batch_size="auto")(\
			# 			delayed(identify_sequence)(file, column_names, \
			#  	[r'GLY_\d+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+GLY_\d+LYS_\d+SER_\d+', \
			#  	r'GLY_\d+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+GLY_\d+LYS_\d+THR_\d+']) for file in files)
			# identify_sequence(files,column_names, \
			# 	[r'GLY_\d+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+GLY_\d+LYS_\d+SER_\d+', \
			# 	r'GLY_\d+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+[A-Z]{3}_[\d]+GLY_\d+LYS_\d+THR_\d+'])
			for file in files:
				matched_list.append(identify_pattern_by_normal(file,pattern, 4))
			pd.DataFrame(matched_list, columns = ['file','pattern_matched'])\
				.to_excel(writer_pattern, sheet_name = "_".join(pattern))

	end = time.time()
	print("Key Search Completed. Total time taken: {} mins".format((end - start)/60))	

	

	
        


