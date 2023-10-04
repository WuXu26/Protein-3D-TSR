
import numpy as np
import pandas as pd
import scipy,csv
import scipy.cluster.hierarchy as hac
from scipy.cluster.hierarchy import ward, average,dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from sklearn.manifold import MDS
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
import seaborn as sns
sns.set_palette("Greens")
from mpl_toolkits.mplot3d import Axes3D

class Dendograming:

	def __init__(self,**kwargs):

		self.setting = kwargs["setting"]
		self.folder = kwargs["outFolder"]
		samples_file ={}
		writer = pd.ExcelWriter(self.folder+'//similarity_values.xlsx', engine='xlsxwriter')
		for name in kwargs["samplesFile"]['sampleClass'].values:
			samples_file[name.split('-')[0]] = name.split('-')[1]

		self.x = pd.read_csv(self.folder+"//normal_jaccard_similarity"+self.setting+".csv",header = 0,index_col=0)
		print self.x.shape
		cols = pd.isnull(self.x).any(1).nonzero()[0]
		
		self.x = self.x.dropna()
		self.x = self.x.drop(self.x.columns[cols],axis=1)
		self.x.columns = list(self.x)
		self.x.index = list(self.x)
		self.x = self.x[sorted(list(self.x))]
		self.x = self.x.sort_index()
		self.similarity_x = (1-self.x) *100
		self.similarity_x.to_excel(writer,sheet_name='Normal Jaccard')
		#self.dist_normal = self.x.values
		
		self.x2 =  pd.read_csv(self.folder+"//generalised_jaccard_similarity"+self.setting+".csv",header = 0,index_col=0)
		self.x2 = self.x2.dropna()
		self.x2 = self.x2.drop(self.x2.columns[cols],axis=1)
		self.x2.columns = list(self.x2)
		self.x2.index = list(self.x2)
		self.x2 = self.x2[sorted(list(self.x2))]
		self.x2 = self.x2.sort_index()
		self.similarity_x2 = (1-self.x2)*100
		self.similarity_x2.to_excel(writer,sheet_name='Generalised Jaccard')
		#self.dist_gen = self.x2.values

		self.x3 =  pd.read_csv(self.folder+"//wu_jaccard_similarity"+self.setting+".csv",header = 0,index_col=0)
		self.x3 = self.x3.dropna()
		self.x3 = self.x3.drop(self.x3.columns[cols],axis=1)
		self.x3.columns = list(self.x3)
		self.x3.index = list(self.x3)
		self.x3 = self.x3[sorted(list(self.x3))]
		self.x3 = self.x3.sort_index()
		self.similarity_x3 = (1-self.x3)*100
		self.similarity_x3.to_excel(writer,sheet_name='Wu')		
		#self.dist_wu = self.x3.values

		self.x4 =  pd.read_csv(self.folder+"//sarika_jaccard1_similarity"+self.setting+".csv",header = 0,index_col=0)
		#
		self.x4 = self.x4.dropna()
		self.x4 = self.x4.drop(self.x4.columns[cols],axis=1)
		self.x4.columns = list(self.x4)
		self.x4.index = list(self.x4)
		self.x4 = self.x4[sorted(list(self.x4))]
		self.x4 = self.x4.sort_index()
		self.similarity_x4 = (1-self.x4)*100
		self.similarity_x4.to_excel(writer,sheet_name='Sarika')	
		#self.similarity_x4.to_csv(self.folder+"//corrected_sarika_jaccard_similarity"+self.setting+".csv")
		#self.dist_sarika1 = self.x4.values

		# self.x1 =  pd.read_csv(self.folder+"//cosine_similarity"+self.setting+".csv",header = 0,index_col=0)
		# self.x1 = self.x1.drop(self.x1.columns[cols],axis=1)
		# self.x1 = self.x1.dropna()
		# self.x1.columns = list(self.x1)
		# self.x1.index = list(self.x1)
		# self.x1 = self.x1[sorted(list(self.x1))]
		# self.x1 = self.x1.sort_index()
		# self.x1.to_excel(writer,sheet_name='Cosine')		
		# #self.x1.to_csv(self.folder+"//corrected_cosine_jaccard_similarity"+self.setting+".csv")
		# self.dist_cosine = self.x1.values

		# self.X = open(kwargs["outFolder"]+'//localFeatureVect'+self.setting +'.csv','r')
		# lines = self.X.readlines()
		# l = []
		# self.X.close()
		# for line in lines:
		# 	i=list(line.split(';')[1].split(','))
		# 	l.append(np.asarray(i[:-1]).astype(np.float))
		# self.X = np.asarray(l)

		self.fileList = sorted(list(self.x))
		self.color_palatte = kwargs["color_palatte"]
		
	def get_hierarchical_clustering_jaccard(self,type):
		#corr_condensed = None
		if type == 'normal':
			corr_condensed = hac.distance.squareform(self.dist_normal)
			#linkage_matrix = average(self.dist_normal) 
		elif type == 'generalised':
			corr_condensed = hac.distance.squareform(self.dist_gen)
			#linkage_matrix = average(self.dist_gen) 
		elif type == 'wu':
			corr_condensed = hac.distance.squareform(self.dist_wu)
			#linkage_matrix = average(self.dist_wu)
		elif type == 'sarika':
			corr_condensed = hac.distance.squareform(self.dist_sarika1)
			#linkage_matrix = average(self.dist_sarika1)
		linkage_matrix = hac.linkage(corr_condensed, method='average')
		
		#Checking the Cophenetic Correclation Coefficient to compare 
		#actual pairwise distances with HClustering implied distance.
		# c, coph_dists = cophenet(linkage_matrix, pdist(self.X))
		# print('Cophenetic correlation: ', c) #This gave 0.336
		# print(hac.fcluster(linkage_matrix,0))
		fig, ax = plt.subplots(figsize=(15, 20)) # set size

		ax = dendrogram(linkage_matrix, orientation="left", labels=self.fileList) 
		plt.tick_params(\
		    axis= 'x',          # changes apply to the x-axis
		    which='both',      # both major and minor ticks are affected
		    labelsize=11,
		    )
		plt.tick_params(\
		    axis= 'y',          # changes apply to the x-axis
		    which='both',      # both major and minor ticks are affected
		    #labelsize=11,
		    )
		plt.title("{} Jaccard Dendogram for {}".format( type,self.setting) )
		plt.savefig('{}/dendo_{}_jaccard{}.png'.format(self.folder,type,self.setting))
		#plt.show()
		plt.gcf().clear()
		
	def get_hierarchical_clustering_cosine(self):
		linkage_matrix = average(self.dist_cosine)  
		#Checking the Cophenetic Correclation Coefficient to compare 
		#actual pairwise distances with HClustering implied distance.
		# c, coph_dists = cophenet(linkage_matrix, pdist(self.X))
		# print('Cophenetic correlation: ', c) #This gave 0.336
		# print(hac.fcluster(linkage_matrix,0))
		fig, ax = plt.subplots(figsize=(15, 20)) # set size

		ax = dendrogram(linkage_matrix, orientation="left", labels=self.fileList) 
		plt.tick_params(\
		    axis= 'x',          # changes apply to the x-axis
		    which='both',      # both major and minor ticks are affected
		    )
		#plt.axvline(x=1.257, c='k')
		#plt.tight_layout() #show plot with tight layout
		plt.title("Cosine Dendogram for " + self.setting)
		plt.savefig(self.folder +'/dendo_cosine'+self.setting+'.png', dpi=200)
		#plt.show()

	def get_heatmap(self, type):
		plt.gcf().clear()
		#cmap = sns.cubehelix_palette(light=8, as_cmap=True)
		#cmap = sns.color_palette("Set1", n_colors=8, desat=.5)
		#cmap = 'summer_r'
		cmap = self.color_palatte
		
		ax = None
		if type == 'normal':
			#ax = sns.heatmap(pd.DataFrame(self.similarity_x), square = True) 
			#linkage = hc.linkage(sp.distance.squareform(self.dist_normal), method='average')
			#ax = sns.clustermap(self.dist_normal, row_linkage=linkage, col_linkage=linkage)
			
			ax = sns.clustermap(pd.DataFrame(self.similarity_x), metric="correlation", cmap=cmap)
			# plt.show()
			# plt.title("{} Jaccard Dendogram for {}".format( type,self.setting) )
			
		if type == 'generalised':
			linkage = hc.linkage(sp.distance.squareform(100 - pd.DataFrame(self.similarity_x2)), method='average')
			ax = sns.clustermap(100 - pd.DataFrame(self.similarity_x2), row_linkage=linkage, col_linkage=linkage)
			#ax = sns.clustermap(pd.DataFrame(self.similarity_x2), metric="correlation", cmap=cmap) #Plot the correlation as heat map
		if type == 'wu':
			ax = sns.clustermap(pd.DataFrame(self.similarity_x3), metric="correlation", cmap=cmap) 
		if type == 'sarika':
			ax =sns.clustermap(pd.DataFrame(self.similarity_x4), metric="correlation", cmap=cmap) 
		if type == 'cosine':
			ax = sns.clustermap(pd.DataFrame(self.dist_cosine), metric="correlation", cmap=cmap) 
		ax.savefig('{}/clustermap_{}_jaccard{}.png'.format(self.folder,type,self.setting))

		#fig = ax.get_figure()
		#fig.savefig('{}/heatmap_{}_{}.png'.format(self.folder,type,self.setting))
		#plt.show()
		plt.gcf().clear()
	def get_heatmap_old(self, type):
		plt.gcf().clear()
		#cmap = 'YlGnBu'
		cmap = self.color_palatte
		ax = None
		if type == 'normal':
			print pd.DataFrame(self.similarity_x)
			ax = sns.heatmap(pd.DataFrame(self.similarity_x), cmap=cmap, square = True, cbar_kws={ "shrink": .8})
		if type == 'generalised':
			ax = sns.heatmap(pd.DataFrame(self.similarity_x2), cmap=cmap, square = True, cbar_kws={ "shrink": .8}) #Plot the correlation as heat map
		if type == 'wu':
			ax = sns.heatmap(pd.DataFrame(self.similarity_x3), cmap=cmap, square = True, cbar_kws={ "shrink": .8})
		if type == 'sarika':
			ax =sns.heatmap(pd.DataFrame(self.similarity_x4), cmap=cmap, square = True, cbar_kws={ "shrink": .8})
		if type == 'cosine':
			ax = sns.heatmap(pd.DataFrame(self.dist_cosine), cmap=cmap, square = True, cbar_kws={ "shrink": .8})

		fig = ax.get_figure()
		fig.savefig('{}/heatmap_{}_{}.png'.format(self.folder,type,self.setting))
		#plt.show()
		plt.gcf().clear()
		
	def get_kmeans_clustering(self,no_of_clusters):
		num_clusters = no_of_clusters
		km = KMeans(n_clusters=num_clusters)
		km.fit(self.X)
		clusters = km.labels_.tolist()
		frame = pd.DataFrame({'protein': self.fileList, 'cluster':clusters}, index = [clusters] , columns = ['protein', 'cluster'])
		#print(frame.sort_values(by = 'cluster'))
		return km,frame  
	  	
	def cluster_analysis(self,frame,no_of_clusters,xs, ys):
		clusters = frame['cluster'].values#km.labels_.tolist()
		colors = ['#1b9e77',  '#d95f02',  '#7570b3', '#e7298a']# , '#66a61e', '#00FFFF',  '#000080',  '#00FF00',  '#FFFF00',  '#808080','#DAF7A6','#581845','#EBDEF0','#F9EBEA','#EAFAF1']
		
		clus =['Cluster 0','Cluster 1','Cluster 2','Cluster 3','Cluster 4']#,'Cluster 5','Cluster 6','Cluster 7','Cluster 8','Cluster 9','Cluster 10','Cluster 11','Cluster 12','Cluster 13','Cluster 14']
		cluster_names = dict(zip(sorted(frame['cluster'].unique()),clus))
		cluster_colors = dict(zip(sorted(frame['cluster'].unique()),colors))
		
		df = pd.DataFrame(dict(x=xs, y=ys, label=clusters, title=frame['protein'].values)) 

		#group by cluster
		groups = df.groupby('label')
		fig, ax = plt.subplots(figsize=(17, 9)) # set size
		ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling

		#iterate through groups to layer the plot
		#note that I use the cluster_name and cluster_color dicts with the 'name' lookup to return the appropriate color/label
		for name, group in groups:
		    ax.plot(group.x, group.y, marker='o', linestyle='', ms=12, 
		            label=cluster_names[name], color=cluster_colors[name], 
		            mec='none')
		    ax.set_aspect('auto')
		    ax.tick_params(\
		        axis= 'x',          # changes apply to the x-axis
		        which='both',      # both major and minor ticks are affected
		        bottom='off',      # ticks along the bottom edge are off
		        top='off',         # ticks along the top edge are off
		        labelbottom='off')
		    ax.tick_params(\
		        axis= 'y',         # changes apply to the y-axis
		        which='both',      # both major and minor ticks are affected
		        left='off',      # ticks along the bottom edge are off
		        top='off',         # ticks along the top edge are off
		        labelleft='off')
		    
		ax.legend(numpoints=1)  #show legend with only 1 point

		#add label in x,y position with the label as the film title
		for i in range(len(df)):
		    ax.text(df.ix[i]['x'], df.ix[i]['y'], df.ix[i]['title'], size=8)  

	    
		#plt.savefig(input_entity_files +'Kmeans_clustering'+setting+setting2+'.png', dpi=200)
		plt.show() #show the plot

	def get_all_figures(self,type):
		plt.gcf().clear()
		cmap = self.color_palatte
		ax = None
		if type == 'normal':
			df = pd.DataFrame(self.similarity_x)
			ax_heatmap = sns.heatmap(df,  square = True, cbar_kws={ "shrink": 0.5})
			linkage = hc.linkage(sp.distance.squareform(100 - df), method='average')
			ax_clustermap = sns.clustermap(100 - df, row_linkage=linkage, col_linkage=linkage)
			#ax_clustermap = sns.clustermap(df, metric="correlation", cmap=cmap)
			corr_condensed = hac.distance.squareform(self.x.values)
		if type == 'generalised':
			ax_heatmap = sns.heatmap(pd.DataFrame(self.similarity_x2),  square = True, cbar_kws={ "shrink": 0.5}) #Plot the correlation as heat map

			linkage = hc.linkage(sp.distance.squareform(100 - pd.DataFrame(self.similarity_x2)), method='average')
			ax_clustermap = sns.clustermap(100 - pd.DataFrame(self.similarity_x2), row_linkage=linkage, col_linkage=linkage)
			corr_condensed = hac.distance.squareform(self.x2.values)
		if type == 'wu':
			df = pd.DataFrame(self.similarity_x3)
			ax_heatmap = sns.heatmap(df,  square = True, cbar_kws={ "shrink": 0.5})
			linkage = hc.linkage(sp.distance.squareform(100 - df), method='average')
			ax_clustermap = sns.clustermap(100 - df, row_linkage=linkage, col_linkage=linkage)
			#ax_clustermap = sns.clustermap(df, metric="correlation", cmap=cmap) 
			corr_condensed = hac.distance.squareform(self.x3.values)
		if type == 'sarika':
			df = pd.DataFrame(self.similarity_x4)
			ax_heatmap =sns.heatmap(df,  square = True, cbar_kws={ "shrink": 0.5})
			linkage = hc.linkage(sp.distance.squareform(100 - df), method='average')
			ax_clustermap = sns.clustermap(100 - df, row_linkage=linkage, col_linkage=linkage)
			#ax_clustermap =sns.clustermap(df, metric="correlation", cmap=cmap)
			corr_condensed = hac.distance.squareform(self.x4.values)
		# if type == 'cosine':
		# 	ax_heatmap = sns.heatmap(pd.DataFrame(self.dist_cosine), cmap=cmap, square = True, cbar_kws={ "shrink": .8})
		# 	ax_clustermap = sns.clustermap(pd.DataFrame(self.dist_cosine), metric="correlation", cmap=cmap) 

		fig = ax_heatmap.get_figure()
		fig.savefig('{}/heatmap_{}_{}.png'.format(self.folder,type,self.setting))

		ax_clustermap.savefig('{}/clustermap_{}_jaccard{}.png'.format(self.folder,type,self.setting))

		linkage_matrix = hac.linkage(corr_condensed, method='average')
		fig, ax = plt.subplots(figsize=(15, 20)) # set size

		ax = dendrogram(linkage_matrix, orientation="left", labels=self.fileList) 
		plt.tick_params(\
		    axis= 'x',          # changes apply to the x-axis
		    which='both',      # both major and minor ticks are affected
		    labelsize=11,
		    )
		plt.tick_params(\
		    axis= 'y',          # changes apply to the x-axis
		    which='both',      # both major and minor ticks are affected
		    #labelsize=11,
		    )
		plt.title("{} Jaccard Dendogram for {}".format( type,self.setting) )
		plt.savefig('{}/dendo_{}_jaccard{}.png'.format(self.folder,type,self.setting))

		plt.gcf().clear()

	def get_dendros_all(self):

		# All Maps
		# Comment any of the below if you have a Memory Error (Sarika)
		self.get_all_figures('normal')
		self.get_all_figures('generalised')
		self.get_all_figures('wu')
		self.get_all_figures('sarika')

		
