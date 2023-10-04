import os
import argparse
from os.path import expanduser
import pandas as pd
import collections
import glob
import time
import ntpath

parser = argparse.ArgumentParser(description='Hierarchical Classification.')
parser.add_argument('--path', '-path', metavar='path', \
    default=os.path.join(expanduser('~'),'Research', 'Protien_Database', \
        'extracted_new_samples', 'testing'), \
    help='Directory of input sample and other files.')
parser.add_argument('--is_header', '-is_header', \
    action='store_true', default=True, \
    help='Enable this option if column header is present on csv files.')
parser.add_argument('--is_append_file_name', '-is_append_file_name', \
    action='store_true', default=True, \
    help='Enable this option if file name is needed in the combined csv.')
parser.add_argument('--column_names', '-column_names', metavar='column_names', \
	default='key,aa0,pos0,aa1,pos1,aa2,pos2,classT,theta,classL,distance,x0,y0,z0,x1,y1,z1,x2,y2,z2',\
	help='Names of the columns of the csv file if they donot have header.')

def merge_csv_without_file_name(all_files):
	#combine all files in the list
	combined = pd.concat([pd.read_csv(f) for f in all_files ])
	return combined

def merge_csv_with_filename(all_files, is_header, column_names):
	all_df = pd.DataFrame(columns = column_names)
	for f in all_files:
		if is_header:
			df = pd.read_csv(f, header = 0)
		else:
			df = pd.read_csv(f, columns = column_names)
		df['file'] = ntpath.basename(f).split('.')[0]
		all_df = pd.concat([all_df, df])
	print(all_df)
	return all_df

if __name__ == '__main__':
	start = time.time()
	args = parser.parse_args()

	outFolder = args.path
	print("Looking for files at: {}".format(outFolder))

	all_files = glob.glob(outFolder+'//*.csv')
	print("Files merging: ")
	print(all_files)	

	column_names = args.column_names.split(',')
	if args.is_append_file_name:		
		combined_csv = merge_csv_with_filename(all_files, args.is_header, column_names)
	else:
		combined_csv = merge_csv_without_file_name(all_files)

	#export to csv
	combined_csv.to_csv( os.path.join(outFolder,"combined_csv.csv"), index=False, encoding='utf-8-sig')

	
	print("Done merging all CSV files")
	print("Time taken for merging(mins):{}".format((time.time()-start)/60))
		
	
