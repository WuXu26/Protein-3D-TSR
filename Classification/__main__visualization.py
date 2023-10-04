
import os
import glob
import ntpath
import socket
import argparse
import time
import pandas as pd

from lib import KeyGeneration, FeatureSelection, Vectorization, JaccardCoefficient, SimilarityScore, Dendograming

parser = argparse.ArgumentParser(description='Parallel Key Generation.')
####Change here for sample_name Sarika - change name in "default = "
parser.add_argument('--sample_name', '-sample', metavar='sample_name', \
	default='sample_hierarchy', help='Name of the sample on which this script should be run.')
parser.add_argument('--thetaBounds', '-theta', metavar='thetaBounds', \
        default='0,12.11,17.32,21.53,25.21,28.54,31.64,34.55,37.34,40.03,42.64,45.17,47.64,50.05,52.43,54.77,57.08,59.38,61.64,63.87,66.09,68.30,70.5,72.69,79.2,81.36,83.51,85.67,87.8,90', \
	help='Bin Boundaries for Theta.')
parser.add_argument('--distBounds', '-dist', metavar='distBounds', 
    default='3.83, 7.00, 9.00, 11.00, 14.00, 17.99, 21.25, 23.19, 24.8, 26.26,27.72, 28.9, 30.36, 31.62, 32.76, 33.84, 35.13, 36.26,37.62,38.73, 40.12,41.8, 43.41, 45.55, 47.46, 49.69, 52.65, 55.81, 60.2, 64.63, 70.04, 76.15,83.26, 132.45', \
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
parser.add_argument('--colors', '-colors', metavar='colors', default=14, \
	help='Select a color palatte number from the list in the commented section below')
# colors = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn',\
# 		 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu',\
# 		 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges',\
# 		  'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', \
# 		  'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', \
# 		  'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', \
# 		  'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', \
# 		  'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', \
# 		  'Spectral_r', 'Vega10', 'Vega10_r', 'Vega20', 'Vega20_r', 'Vega20b', 'Vega20b_r', \
# 		  'Vega20c', 'Vega20c_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', \
# 		  'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', \
# 		  'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cool', 'cool_r', \
# 		  'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', \
# 		  'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', \
# 		  'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', \
# 		  'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', \
# 		  'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r', 'inferno', 'inferno_r', 'jet', 'jet_r', \
# 		  'magma', 'magma_r', 'mako', 'mako_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', \
# 		  'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', \
# 		  'rocket_r', 'seismic', 'seismic_r', 'spectral', 'spectral_r', 'spring', 'spring_r', 'summer', \
# 		  'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', \
# 		  'terrain', 'terrain_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r', 'winter', 'winter_r']




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

##	# # Feature Selection	
##	print('--------------------------Start FeatureSelection-------------------------------')
##	features = FeatureSelection(outFolder = outFolder, setting = setting, numOfGap = kwargs['numGap'],mad = kwargs['mad'], \
##			filesList = filesList, keyCombine = kwargs['keyCombine'],featureSelection = kwargs['featureSelection'] )
##	features.feature_selection()
##	print('--------------------------End FeatureSelection----------------------------------')
##
##	# # Vectorization
##	print('--------------------------Start Vectorization-----------------------------------')
##	if kwargs['keyCombine'] == 0:
##		changedFiles=glob.glob(kwargs["outFolder"]+'//*.keys_'+kwargs['outFolderName'])
##	else:
##		changedFiles=glob.glob(kwargs["outFolder"]+'//*.keys_keycombine'+str(kwargs['keyCombine']))
##	changedList = sorted([ntpath.basename(file) for file in changedFiles])
##	
##	vectors = Vectorization(outFolder = outFolder, setting = setting, filesList = changedList)
##	vectors.vectorize()
##	print('--------------------------End Vectorization-------------------------------------')
##
##
##	# JaccardCoefficient
##	print('--------------------------Start Jaccard-----------------------------------------')
##	jaccard = JaccardCoefficient(outFolder = outFolder, setting = setting, filesList = filesList, normalize=kwargs['normalize'],sample_dict = df_dict)
##	jaccard.calculate_jaccard()
##	print('--------------------------End Jaccard-------------------------------------------')

	# Dendogram
	print('--------------------------Start Clustering--------------------------------------')
	dendo = Dendograming(samplesFile = df, outFolder = outFolder, setting = setting, filesClass = fileClass , color_palatte = kwargs['color_palatte'])
	dendo.get_dendros_all()
	print('--------------------------End Clustering----------------------------------------')


if __name__ == '__main__':
	"""Executable code starts here."""
	args = parser.parse_args()
	if socket.gethostname().startswith('qb'):
		path = '/work/wxx6941/TSR/Protein_Database/extracted_samples/testing/'
	else:
		path = '/work/wxx6941/TSR/Protein_Database/extracted_samples/testing/'
	subFolder= args.sample_name + '/'
	sampleDetailsFile = path+subFolder+'sample_details.csv'
	thetaBounds = list(map(float, args.thetaBounds.split(',')))
	distBounds =list(map(float, args.distBounds.split(',')))

	outFolderName = 'theta'+str(len(thetaBounds) - 1)+'_dist'+str(len(distBounds) + 1)
	outFolder = path +subFolder+outFolderName
	amino_acid_code=None
	#amino_acid_code=open(path+"aminoAcidCode_lexicographic _new.txt","r") 

	start_time=time.time()
	colors = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn',\
		 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu',\
		 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges',\
		  'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', \
		  'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', \
		  'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', \
		  'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', \
		  'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', \
		  'Spectral_r', 'Vega10', 'Vega10_r', 'Vega20', 'Vega20_r', 'Vega20b', 'Vega20b_r', \
		  'Vega20c', 'Vega20c_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', \
		  'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', \
		  'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cool', 'cool_r', \
		  'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', \
		  'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', \
		  'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', \
		  'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', \
		  'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r', 'inferno', 'inferno_r', 'jet', 'jet_r', \
		  'magma', 'magma_r', 'mako', 'mako_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', \
		  'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', \
		  'rocket_r', 'seismic', 'seismic_r', 'spectral', 'spectral_r', 'spring', 'spring_r', 'summer', \
		  'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', \
		  'terrain', 'terrain_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r', 'winter', 'winter_r']



	run_classification_pipeline(path = path, subFolder= subFolder, thetaBounds = thetaBounds, outFolderName= outFolderName, sample_name=args.sample_name,\
			distBounds=distBounds, outFolder=outFolder, filesType=args.filesType, amino_acid_code=amino_acid_code, featureSelection = args.featureSelection,\
			normalize = args.normalize, numGap = args.numGap, mad = args.mad, keyCombine = args.keyCombine, numOfLabels= args.numOfLabels, sampleDetailsFile = sampleDetailsFile , color_palatte = colors[int(args.colors)])

	end_time=time.time()
	total_time=((end_time)-(start_time))
	print("Classification done. Total Time taken(secs): {}".format(total_time))
