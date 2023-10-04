
import os
import glob
import ntpath
import socket
import argparse
import time
import pandas as pd

from lib import KeyGeneration, FeatureSelection, Vectorization, JaccardCoefficient, SimilarityScore, Dendograming

parser = argparse.ArgumentParser(description='Parallel Key Generation.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', \
	default='sample_protease', help='Name of the sample on which this script should be run.')
parser.add_argument('--thetaBounds', '-theta', metavar='thetaBounds', \
	default='0, 14.1, 20.31, 25.49, 29.99, 34.1, 38, 41.7, 45.3, 48.61, 51.91, 55.19, 58.29, 61.3, 64.39, 67.3, 70.3, 73.11, 81.69, 84.49, 87.29, 90', \
	help='Bin Boundaries for Theta.')
parser.add_argument('--distBounds', '-dist', metavar='distBounds', 
	default='2.81, 7.00, 9.00, 11.00, 14.00, 17.4, 24.16, 30.19, 36.37, 44.78, 175.52', \
	help='Bin Boundaries for maxDist.')
parser.add_argument('--filesType', '-filesType', metavar='filesType', \
	default='*.pdb', help='Type of the downloaded protein file.')  # "*.ent"
parser.add_argument('--featureSelection', '-featureSelection', metavar='featureSelection', \
	default=False, help='Argument to tell if Feature Selection is required.')
parser.add_argument('--normalize', '-normalize', metavar='normalize', default=False, \
	help='Argument to tell if Normalization is required.')
parser.add_argument('--numGap', '-numGap', metavar='numGap', default=9, \
	help='Feature Selection related parameter.')
parser.add_argument('--mad', '-mad', metavar='mad', default=0, \
	help='Feature Selection related parameter.')
parser.add_argument('--keyCombine', '-keyCombine', metavar='keyCombine', default=0, \
	help='Argument to deal with higher level grouping of keys.')
parser.add_argument('--numOfLabels', '-numOfLabels', metavar='numOfLabels', default=20, \
	help='Set it to 12 without Amino Acid grouping.')
parser.add_argument('--normalJaccard', '-normalJaccard', metavar='normalJaccard', default=True, \
	help='Set this to False if you donot need Normal Jaccard similarity calculations.')
parser.add_argument('--generalisedJaccard', '-generalisedJaccard', metavar='generalisedJaccard', default=True, \
	help='Set this to False if you donot need Generalized Jaccard similarity calculations.')
parser.add_argument('--wuJaccard', '-wuJaccard', metavar='wuJaccard', default=True, \
	help='Set this to False if you donot need Wu Generalised Jaccard similarity calculations.')
parser.add_argument('--sarikaJaccard', '-sarikaJaccard', metavar='sarikaJaccard', default=True, \
	help='Set this to False if you donot need Sarika Jaccard similarity calculations.')


def run_classification_pipeline(**kwargs):
	
	setting = '_' + kwargs['outFolderName'] + '_gap' + str(kwargs['numGap']) + \
			  '_mad' + str(kwargs['mad']) + '_keyCombine' + str(kwargs['keyCombine']) \
			  if kwargs['featureSelection'] else '_' + outFolderName + '_NoFeatureSelection' + '_keyCombine' + str(kwargs['keyCombine'])
	print(setting)
	print("Working on sample {} theta bins: {} length bins: {}". \
			  format(kwargs['sample_name'], str(len(kwargs['thetaBounds']) - 1), str(len(kwargs['distBounds']) + 1)))
	# Key Generation
	# files=glob.glob(kwargs["path"]+kwargs["subFolder"]+kwargs["filesType"])
	# keys = KeyGeneration(path = kwargs["path"],subFolder = kwargs["subFolder"],fileType = kwargs["filesType"] ,aminoAcidCode = AMINO_ACID_CODE, thetaBounds = THETA_BOUNDS, distBounds = DIST_BOUNDS, numOfLabels = NUM_LABELS)
	# Parallel(n_jobs=cpu_count() - 1, verbose=10, backend="multiprocessing", batch_size="auto")(delayed(keys.processFiles)(fileName) for fileName in files)
	# for fileName in files:
	#	keys.processFiles(fileName)
	# print("Key Generation Complete.")

	os.chdir(outFolder)
	files=glob.glob(kwargs["outFolder"]+'//*.keys_'+kwargs['outFolderName'])	
	filesList = sorted([ntpath.basename(file) for file in files])	

	df = pd.read_csv(kwargs["sampleDetailsFile"])
	df['protein'] = df['protein'].apply(lambda x: x.upper())
	df['sampleClass'] = df['group'] +'-'+df['protein']
	df['sampleClass'] = map(lambda x: x.upper(), df['sampleClass'])
	fileClass = df['sampleClass'].values
	df_dict = dict(zip(df.protein,df.group))

	# # Feature Selection	
	print('--------------------------Start FeatureSelection-------------------------------')
	features = FeatureSelection(outFolder = outFolder, setting = setting, numOfGap = kwargs['numGap'],mad = kwargs['mad'], \
			filesList = filesList, keyCombine = kwargs['keyCombine'],featureSelection = kwargs['featureSelection'] )
	features.feature_selection()
	print('--------------------------End FeatureSelection----------------------------------')

	# # Vectorization
	print('--------------------------Start Vectorization-----------------------------------')
	if kwargs['keyCombine'] == 0:
		changedFiles=glob.glob(kwargs["outFolder"]+'//*.keys_'+kwargs['outFolderName'])
	else:
		changedFiles=glob.glob(kwargs["outFolder"]+'//*.keys_keycombine'+str(kwargs['keyCombine']))
	changedList = sorted([ntpath.basename(file) for file in changedFiles])
	
	vectors = Vectorization(outFolder = outFolder, setting = setting, filesList = changedList)
	vectors.vectorize()
	print('--------------------------End Vectorization-------------------------------------')


	# JaccardCoefficient
	print('--------------------------Start Jaccard-----------------------------------------')
	jaccard = JaccardCoefficient(outFolder = outFolder, setting = setting, filesList = filesList, normalize=kwargs['normalize'],sample_dict = df_dict)
	jaccard.calculate_jaccard()
	print('--------------------------End Jaccard-------------------------------------------')

	# Dendogram
	print('--------------------------Start Clustering--------------------------------------')
	dendo = Dendograming(samplesFile = df, outFolder = outFolder, setting = setting, filesClass = fileClass)
	dendo.get_dendros_all()
	print('--------------------------End Clustering----------------------------------------')


if __name__ == '__main__':
	"""Executable code starts here."""
	args = parser.parse_args()
	if socket.gethostname().startswith('qb'):
		path = '//work//wxx6941//TSR//Protein_Database//'
	else:
		path = '/home/linc/c00219805/Research/Protien_Database/'
	subFolder= '//extracted_samples/testing/'+args.sample_name+'//'
	sampleDetailsFile = path+subFolder+'sample_details.csv'
	thetaBounds = list(map(float, args.thetaBounds.split(',')))
	distBounds =list(map(float, args.distBounds.split(',')))

	outFolderName = 'theta'+str(len(thetaBounds) - 1)+'_dist'+str(len(distBounds) + 1)
	outFolder = path +subFolder+outFolderName
	amino_acid_code=open(path+"aminoAcidCode_lexicographic_new.txt","r") 

	start_time=time.time()
	

	run_classification_pipeline(path = path, subFolder= subFolder, thetaBounds = thetaBounds, outFolderName= outFolderName, sample_name=args.sample_name,\
			distBounds=distBounds, outFolder=outFolder, filesType=args.filesType, amino_acid_code=amino_acid_code, featureSelection = args.featureSelection,\
			normalize = args.normalize, numGap = args.numGap, mad = args.mad, keyCombine = args.keyCombine, numOfLabels= args.numOfLabels, sampleDetailsFile = sampleDetailsFile)

	end_time=time.time()
	total_time=((end_time)-(start_time))
	print("Classification done. Total Time taken(secs): {}".format(total_time))
