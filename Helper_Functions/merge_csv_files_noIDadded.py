import os
import argparse
from os.path import expanduser
import pandas as pd
import collections
import glob
import time

parser = argparse.ArgumentParser(description='Hierarchical Classification.')
parser.add_argument('--path', '-path', metavar='path', \
    default=os.path.join(expanduser('~'),'Research', 'Protien_Database', \
        'extracted_new_samples', 'testing'), \
    help='Directory of input sample and other files.')



if __name__ == '__main__':
	start = time.time()
	args = parser.parse_args()

	outFolder = args.path
	print("Looking for files at: {}".format(outFolder))

	all_files = glob.glob(outFolder+'//*_key_search_common_keys.csv')
	print("Files merging: ")
	print(all_files)	

	#combine all files in the list
	combined_csv = pd.concat([pd.read_csv(f) for f in all_files ])
	#export to csv
	combined_csv.to_csv( os.path.join(outFolder,"combined_csv.csv"), index=False, encoding='utf-8-sig')
	
	print("Done merging all CSV files")
	print("Time taken for merging(mins):{}".format((time.time()-start)/60))

		
	
