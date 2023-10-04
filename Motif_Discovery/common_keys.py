import glob, os, csv, ntpath,socket,argparse
import pandas as pd, numpy as np
from random import randint

parser = argparse.ArgumentParser(description='Generate similar, identical, disimilar and all key details files and do analysis on them.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', default='sample_protease', help='Name of the sample on which this script should be run.')
parser.add_argument('--run_preprocessing', '-preprocess', metavar='run_preprocessing', default=False, help='Argument to enable pre- processing steps.')
parser.add_argument('--run_postprocessing', '-postprocess', metavar='run_postprocessing', default=False, help='Argument to enable post- processing steps.')
parser.add_argument('--run_random100', '-random100', metavar='run_random100', default=False, help='Argument to extract random 100 samples of the pre-processing output files.')
parser.add_argument('--thetaBins', '-theta', metavar='thetaBins', default=21, help='Number of bins for Theta.')
parser.add_argument('--distBins', '-dist', metavar='distBins', default=12, help='Number of bins for maxDist.')      


def processing():
	calculationsRequired = []
	file1Details = []
	file2Details =[]
	
	df1_keys_grp = pd.read_table(path + subfolder+bins+'/'+file1+'.keys_original_keycombine2' , delim_whitespace=True)
	df1_keys_grp['fileName']= file1
	
	df2_keys_grp = pd.read_table(path + subfolder+bins+'/'+file2+'.keys_original_keycombine2' , delim_whitespace=True)
	df2_keys_grp['fileName']= file2
		
	#print len(df1_keys_grp),len(df1),len(df2_keys_grp),len(df2)
	calculationsRequired.append('No. of original keys')
	file1Details.append(df1_keys_grp['freq'].sum())
	file2Details.append(df2_keys_grp['freq'].sum())

	#Obtaing all the common new keys in both the files
	common_newkeys = list( set(list(set(df1_keys_grp['grp'].values) & set(df2_keys_grp['grp'].values))))
	all_newkeys = list( set(list(set(df1_keys_grp['grp'].values) | set(df2_keys_grp['grp'].values))))
	outside_newkeys = list(set(list(set(all_newkeys) - set(common_newkeys))))
	df_appened_dissimilar = df1_keys_grp[df1_keys_grp['grp'].isin(outside_newkeys)].append(df2_keys_grp[df2_keys_grp['grp'].isin(outside_newkeys)], ignore_index=True).sort_values(by = ['key'])
	calculationsRequired.append('No. of dissimilar keys')
	file1Details.append(df1_keys_grp[df1_keys_grp['grp'].isin(outside_newkeys)]['freq'].sum())
	file2Details.append(df2_keys_grp[df2_keys_grp['grp'].isin(outside_newkeys)]['freq'].sum())
	
	print('all- ', len(all_newkeys), 'similar- ', len(common_newkeys), 'dissimilar- ', len(outside_newkeys))
	df11= df1_keys_grp[df1_keys_grp['grp'].isin(common_newkeys)]
	df22 = df2_keys_grp[df2_keys_grp['grp'].isin(common_newkeys)]

	all_keys = list( set(list(set(df1_keys_grp['key'].values) | set(df2_keys_grp['key'].values))))
	identical_keys = list( set(list(set(df11['key'].values) & set(df22['key'].values))))
	
	df_similar = df11.append(df22, ignore_index=True).sort_values(by = ['key'])
	calculationsRequired.append('No. of similar original keys')
	file1Details.append(df11['freq'].sum())
	file2Details.append(df22['freq'].sum())

	df_appened_identical = df11[df11['key'].isin(identical_keys)].append(df22[df22['key'].isin(identical_keys)], ignore_index=True).sort_values(by = ['key'])
	calculationsRequired.append('No. of identical keys')
	file1Details.append(df11[df11['key'].isin(identical_keys)]['freq'].sum())
	file2Details.append(df22[df22['key'].isin(identical_keys)]['freq'].sum())

	print('Started writing to excel..')

	writer = pd.ExcelWriter(path+subfolder+bins+'2-protein-plus_minus_2_analysis.xlsx', engine='xlsxwriter')

	df_summary= pd.DataFrame({'values':calculationsRequired,file1:file1Details,file2:file2Details})
	print(df_summary)
	df_summary.to_excel(writer,sheet_name='Summary')
	df_similar.groupby(['key','grp','fileName'])['freq'].agg('sum').reset_index(name='sum').to_excel(writer, sheet_name='similar_keys')
	df_appened_identical.groupby(['key','grp','fileName'])['freq'].agg('sum').reset_index(name='sum').to_excel(writer, sheet_name='identicalKeys')
	df_appened_dissimilar.groupby(['key','grp','fileName'])['freq'].agg('sum').reset_index(name='sum').to_excel(writer, sheet_name='dissimilarKeys')
	df1_keys_grp.to_excel(writer,sheet_name='all_keys_'+file1)
	df2_keys_grp.to_excel(writer,sheet_name='all_keys_'+file2)
	writer.save()

	print('Started writing to CSV..')
	df1_triplets = pd.read_table(path + subfolder+bins+'/'+file1+'.triplets_theta21_dist12' , delim_whitespace=True,names = ('key','aacd0','position0','aacd1','position1','aacd2','position2','classT1','Theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'))
	df1_triplets = df1_triplets[['key','aacd0','position0','aacd1','position1','aacd2','position2','classT1','Theta','classL1','maxDist']]
	df1 = df1_keys_grp.merge(df1_triplets,on = 'key',how='left')

	df2_triplets = pd.read_table(path + subfolder+bins+'/'+file2+'.triplets_theta21_dist12' , delim_whitespace=True,names = ('key','aacd0','position0','aacd1','position1','aacd2','position2','classT1','Theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'))
	df2_triplets = df2_triplets[['key','aacd0','position0','aacd1','position1','aacd2','position2','classT1','Theta','classL1','maxDist']]
	df2 = df2_keys_grp.merge(df2_triplets,on = 'key',how='left')

	#For size purposes separated the two files. If no size issues comment individual file writing and uncomment appending code.
	df1[df1['key'].isin(df_similar['key'])].to_csv(path+subfolder+bins+'/similarKeys_'+file1+'_'+setting +'.csv', sep=',')
	df2[df2['key'].isin(df_similar['key'])].to_csv(path+subfolder+bins+'/similarKeys_'+file2+'_'+setting +'.csv', sep=',')
	#df1[df1['key'].isin(df_similar['key'])].append(df2[df2['key'].isin(df_similar['key'])]).to_csv(path+subfolder+bins+'/similarKeys_'+file1+'_'+file2+'_'+setting +'.csv', sep=',')
	df1[df1['key'].isin(df_appened_identical['key'])].to_csv(path+subfolder+bins+'/identicalKeys_'+file1+'_'+setting +'.csv', sep=',')
	df2[df2['key'].isin(df_appened_identical['key'])].to_csv(path+subfolder+bins+'/identicalKeys_'+file2+'_'+setting +'.csv', sep=',')
	#df1[df1['key'].isin(df_appened_identical['key'])].append(df2[df2['key'].isin(df_appened_identical['key'])]).to_csv(path+subfolder+bins+'/identicalKeys_'+file1+'_'+file2+'_'+setting +'.csv', sep=',')
	df1[df1['key'].isin(df_appened_dissimilar['key'])].to_csv(path+subfolder+bins+'/dissimilarKeys_'+file1+'_'+setting +'.csv', sep=',')
	df2[df2['key'].isin(df_appened_dissimilar['key'])].to_csv(path+subfolder+bins+'/dissimilarKeys_'+file2+'_'+setting +'.csv', sep=',')
	#df1[df1['key'].isin(df_appened_dissimilar['key'])].append(df2[df2['key'].isin(df_appened_dissimilar['key'])]).to_csv(path+subfolder+bins+'/dissimilarKeys_'+file1+'_'+file2+'_'+setting +'.csv', sep=',')
	#merged.to_csv(path+subfolder+bins+'/dissimilarKeys_with_fileFreqs_'+file1+'_'+file2+'_'+setting +'.csv', sep=',')
	print('Finished writing to CSV..')

	print('Done.')

def generateFrequencyGroups(x):
	x = int(x)	
	group = ''
	if x <10:
		group = '1-10'
	elif 11 <= x < 100:
		group = '11-100'
	elif 101 <= x < 500:
		group = '101-500'
	else:
		group = '>500'
	return group

def calculateAminoAcidDistribution(keysFile1,keysFile2,isFullCount,isGroupFreq):
	'''Calculates the percent of amino acids in different types of keys file.'''

	all_Keys = keysFile1.append(keysFile2)
	groubyColumns = ['aacd','fileName']

	if isGroupFreq:
		all_Keys['freqGroups'] = all_Keys['freq'].apply(generateFrequencyGroups)
		aacd0 = all_Keys[['freq','aacd0','fileName','freqGroups']]
		aacd0.columns= ['freq','aacd','fileName','freqGroups']
		aacd1 = all_Keys[['freq','aacd1','fileName','freqGroups']]
		aacd1.columns= ['freq','aacd','fileName','freqGroups']
		aacd2 = all_Keys[['freq','aacd2','fileName','freqGroups']]
		aacd2.columns= ['freq','aacd','fileName','freqGroups']
		groubyColumns = ['aacd','fileName','freqGroups']

	else:	
		aacd0 = all_Keys[['freq','aacd0','fileName']]
		aacd0.columns= ['freq','aacd','fileName']
		aacd1 = all_Keys[['freq','aacd1','fileName']]
		aacd1.columns= ['freq','aacd','fileName']
		aacd2 = all_Keys[['freq','aacd2','fileName']]
		aacd2.columns= ['freq','aacd','fileName']
		
	
	aminoacid_freq = aacd0.append(aacd1).append(aacd2)

	#If you want to count frequencies only once isFullCount = False.
	if not isFullCount:
		aminoacid_freq.loc[aminoacid_freq['freq']>1,'freq'] = 1

	aminoacid_freq = aminoacid_freq.groupby(groubyColumns)['freq'].agg('sum').reset_index(name='sum')
			

	aminoacid_freq_file1 = aminoacid_freq[aminoacid_freq['fileName'] == file1]	
	aminoacid_freq_file2 = aminoacid_freq[aminoacid_freq['fileName'] == file2]

	#If youwant amino acid distribution per frequency group.
	#Frequency Groups are mentioned in generateFrequencyGroups function.
	if isGroupFreq:
		file1_dict = aminoacid_freq_file1[['sum','freqGroups']].groupby(['freqGroups']).agg('sum').to_dict()
		aminoacid_freq_file1['totals'] = aminoacid_freq_file1['freqGroups'].apply(lambda x: file1_dict['sum'][x])
		aminoacid_freq_file1['percent'] = 100*aminoacid_freq_file1['sum']/aminoacid_freq_file1['totals']

		file2_dict = aminoacid_freq_file2[['sum','freqGroups']].groupby(['freqGroups']).agg('sum').to_dict()
		aminoacid_freq_file2['totals'] = aminoacid_freq_file2['freqGroups'].apply(lambda x: file2_dict['sum'][x])
		aminoacid_freq_file2['percent'] = 100*aminoacid_freq_file2['sum']/aminoacid_freq_file2['totals']
	else:
		aminoacid_freq_file1['percent'] = aminoacid_freq_file1['sum']/aminoacid_freq_file1['sum'].sum()
		aminoacid_freq_file2['percent'] = aminoacid_freq_file2['sum']/aminoacid_freq_file2['sum'].sum()

	#print(aminoacid_freq_file1.append(aminoacid_freq_file2))
	
	return aminoacid_freq_file1.append(aminoacid_freq_file2)

def calculateMeanSD(file,fileType,isTripletFile):
	if not isTripletFile:
		df = pd.read_csv(path+subfolder+bins+'/'+fileType+'_'+file+'_'+setting +'.csv')
	else:
		df = pd.read_table(path + subfolder+bins+'/'+file+'.triplets_theta21_dist12' , delim_whitespace=True,names = ('key','aacd0','position0','aacd1','position1',
		'aacd2','position2','classT1','Theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'))
	return file,fileType,df['Theta'].mean(),df['Theta'].std(),df['maxDist'].mean(),df['maxDist'].std()


def postProcessing():
	"""Calculate Mean and SD of Theta and maxDist. Calculate amino acid distribution of all file categories."""

	#Write to Excel File.###
	writer = pd.ExcelWriter(path+subfolder+bins+'/'+file1+'_'+file2+'_Mean_AA_distribution_analysis.xlsx', engine='xlsxwriter')

	#Calculate Mean and Standand Deviation of Thetas and maxDists
	alltuples = []
	alltuples.append(calculateMeanSD(file1,'similarKeys',False))
	alltuples.append(calculateMeanSD(file2,'similarKeys',False))
	alltuples.append(calculateMeanSD(file1,'identicalKeys',False))
	alltuples.append(calculateMeanSD(file2,'identicalKeys',False))
	alltuples.append(calculateMeanSD(file1,'dissimilarKeys',False))
	alltuples.append(calculateMeanSD(file2,'dissimilarKeys',False))
	alltuples.append(calculateMeanSD(file1,'AllKeys',True))
	alltuples.append(calculateMeanSD(file2,'AllKeys',True))
	pd.DataFrame(alltuples,columns=['fileName','fileType' ,'mean_theta', 'SD_theta','mean_maxDist','SD_maxDist']).to_excel(writer,sheet_name='Mean and SD')
	print('Mean and SD done.')

	#calculate amino acid percent from similar keys files.
	similarKeys_file1 = pd.read_csv(path+subfolder+bins+'/similarKeys_'+file1+'_'+setting +'.csv')[['key','freq','grp','fileName','aacd0','aacd1','aacd2']].drop_duplicates()
	similarKeys_file2 = pd.read_csv(path+subfolder+bins+'/similarKeys_'+file2+'_'+setting +'.csv')[['key','freq','grp','fileName','aacd0','aacd1','aacd2']].drop_duplicates()
	calculateAminoAcidDistribution(similarKeys_file1,similarKeys_file2,True,True).to_excel(writer,sheet_name='Similar_Grouping')
	calculateAminoAcidDistribution(similarKeys_file1,similarKeys_file2,False,True).to_excel(writer,sheet_name='Freq1_Similar_Grouping')
	calculateAminoAcidDistribution(similarKeys_file1,similarKeys_file2,True,False).to_excel(writer,sheet_name='Similar')
	calculateAminoAcidDistribution(similarKeys_file1,similarKeys_file2,False,False).to_excel(writer,sheet_name='Freq1_Similar')
	print('Similar keys Amino Acid distribution done.')

	#calculate amino acid percent from identical keys files.
	identicalKeys_file1 = pd.read_csv(path+subfolder+bins+'/identicalKeys_'+file1+'_'+setting +'.csv')[['key','freq','grp','fileName','aacd0','aacd1','aacd2']].drop_duplicates()
	identicalKeys_file2 = pd.read_csv(path+subfolder+bins+'/identicalKeys_'+file2+'_'+setting +'.csv')[['key','freq','grp','fileName','aacd0','aacd1','aacd2']].drop_duplicates()
	calculateAminoAcidDistribution(identicalKeys_file1,identicalKeys_file2,True,True).to_excel(writer,sheet_name='Identical_Grouping')
	calculateAminoAcidDistribution(identicalKeys_file1,identicalKeys_file2,False,True).to_excel(writer,sheet_name='Freq1_Identical_Grouping')
	calculateAminoAcidDistribution(identicalKeys_file1,identicalKeys_file2,True,False).to_excel(writer,sheet_name='Identical')
	calculateAminoAcidDistribution(identicalKeys_file1,identicalKeys_file2,False,False).to_excel(writer,sheet_name='Freq1_Identical')
	print('Identical keys Amino Acid distribution done.')

	#calculate amino acid percent from dissimilar keys files.
	dissimilarKeys__file1 = pd.read_csv(path+subfolder+bins+'/dissimilarKeys_'+file1+'_'+setting +'.csv')[['key','freq','grp','fileName','aacd0','aacd1','aacd2']].drop_duplicates()
	dissimilarKeys__file2 = pd.read_csv(path+subfolder+bins+'/dissimilarKeys_'+file2+'_'+setting +'.csv')[['key','freq','grp','fileName','aacd0','aacd1','aacd2']].drop_duplicates()
	calculateAminoAcidDistribution(dissimilarKeys__file1,dissimilarKeys__file2,True,True).to_excel(writer,sheet_name='Dissimilar_Grouping')
	calculateAminoAcidDistribution(dissimilarKeys__file1,dissimilarKeys__file2,False,True).to_excel(writer,sheet_name='Freq1_Dissimilar_Grouping')
	calculateAminoAcidDistribution(dissimilarKeys__file1,dissimilarKeys__file2,True,False).to_excel(writer,sheet_name='Dissimilar')
	calculateAminoAcidDistribution(dissimilarKeys__file1,dissimilarKeys__file2,False,False).to_excel(writer,sheet_name='Freq1_Dissimilar')
	print('Dissimilar keys Amino Acid distribution done.')

	#calculate amino acid percent from all keys files.
	df1_triplets = pd.read_table(path + subfolder+bins+'/'+file1+'.triplets_theta21_dist12' , delim_whitespace=True,names = ('key','aacd0','position0','aacd1','position1',
		'aacd2','position2','classT1','Theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'))[['key','aacd0','aacd1','aacd2']].drop_duplicates()
	df1 = pd.read_table(path + subfolder+bins+'/'+file1+'.keys_original_keycombine2' , delim_whitespace=True).merge(df1_triplets,on = 'key',how='left')
	df1['fileName']= file1
	df2_triplets = pd.read_table(path + subfolder+bins+'/'+file2+'.triplets_theta21_dist12' , delim_whitespace=True,names = ('key','aacd0','position0','aacd1','position1',
		'aacd2','position2','classT1','Theta','classL1','maxDist','x0','y0','z0','x1','y1','z1','x2','y2','z2'))[['key','aacd0','aacd1','aacd2']].drop_duplicates()
	df2 = pd.read_table(path + subfolder+bins+'/'+file2+'.keys_original_keycombine2' , delim_whitespace=True).merge(df2_triplets,on = 'key',how='left')
	df2['fileName']= file2
	calculateAminoAcidDistribution(df1,df2,True,True).to_excel(writer,sheet_name='All Keys_Grouping')
	calculateAminoAcidDistribution(df1,df2,False,True).to_excel(writer,sheet_name='Freq1_All Keys_Grouping')
	calculateAminoAcidDistribution(df1,df2,True,False).to_excel(writer,sheet_name='All Keys')
	calculateAminoAcidDistribution(df1,df2,False,False).to_excel(writer,sheet_name='Freq1_All Keys')
	print('All keys Amino Acid distribution done.')	
	
	writer.save()


def randomSelectionKeys(file1,file2):
	file1 = pd.read_csv(file1 )[['fileName','key','aacd0','position0','aacd1','position1','aacd2','position2','classT1','Theta','classL1','maxDist']]
	randomRows = [ randint(1, len(file1)-1) for x in range(100)]
	print(randomRows, len(randomRows))
	alltuples = [file1.loc[rowNo] for rowNo in randomRows]
		
	df_file1= pd.DataFrame(alltuples,columns=['fileName','key','aacd0','position0','aacd1','position1','aacd2','position2','classT1','Theta','classL1','maxDist'])
	print(df_file1)
	identifiedKeys = df_file1['key'].values
	df_file1 = file1[file1['key'].isin(identifiedKeys)]

	file2 = pd.read_csv(file2 )[['fileName','key','aacd0','position0','aacd1','position1','aacd2','position2','classT1','Theta','classL1','maxDist']]
	df_file2 = file2[file2['key'].isin(identifiedKeys)]
	print(df_file2)


	writer = pd.ExcelWriter(path+subfolder+bins+'random100_identical_keys.xlsx', engine='xlsxwriter')
	pd.DataFrame({'Random Keys': identifiedKeys}).to_excel(writer,sheet_name='100 Keys')
	df_file1.append(df_file2).to_excel(writer,sheet_name='Random100')
	writer.save()
	
if __name__ == '__main__':
	"""Executable code starts here."""

	args = parser.parse_args()

	#Add more conditions if you are running on different servers.
	if socket.gethostname() == 'qb2':
	    path = '/work/sarikak/TSR/Protien_Database/'
	else:
	    path = '/home/linc/c00219805/Research/Protien_Database/'
	subfolder= 'extracted_new_samples/testing/'+args.sample_name+'/'
	bins = '/theta{}_dist{}/'.format(args.thetaBins,args.distBins)
	files=glob.glob(path+subfolder+"*.pdb")
	fileNames= [ntpath.basename(i)[:-4] for i in files]
	file1=fileNames[0]
	file2= fileNames[1]

	fileType = "*.keys_keycombine2"  
	fileType_triplets = "*.triplets_theta21_dist12" 
	setting = 'keycombine2'
	files_triplets = glob.glob(path+subfolder+bins+fileType_triplets)
	print(fileNames)

	if args.run_preprocessing:
		print('Started Pre-Processing..')
		processing()
		print('Pre-Processing Done.')
	if args.run_postprocessing:
		print('Started Post-Processing..')
		postProcessing()
		print('Post-Processing Done.')
	#randomSelectionKeys(path+subfolder+bins+'/identicalKeys_'+file1+'_'+setting+'.csv', path+subfolder+bins+'/identicalKeys_'+file2+'_'+setting+'.csv')

