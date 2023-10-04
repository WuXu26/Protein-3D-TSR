import numpy as np
import pandas as pd
import scipy,csv
import scipy.cluster.hierarchy as hac
from scipy.cluster.hierarchy import ward, average,dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances
from sklearn.manifold import MDS
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
import gensim
from gensim import corpora

setting = '_all_features'# '_feature_selected'#'_gap0_nomad'#'_theta16_dist16_gap25_mad0'
#setting2 = '_jaccard_kmeans'
setting2 = '_cosine_kmeans'
fileList=['1HCK-CMGC','1GNG-CMGC','1HOW-CMGC','1JNK-CMGC','1LUF-TK','2SRC-TK','1M14-TK','1IR3-TK','1GJO-TK','1O6L-AGC','1H1W-AGC','1CDK-AGC','1OMW-AGC','1A06-CAMK','1KWP-CAMK','1TKI-CAMK','1JKL-CAMK','1PHK-CAMK','1IA8-CAMK','1O6Y-Other','1MUO-Other','1LP4-Other','1CJA-Atyp','1BO1-Atyp','1NW1-Atyp','1J7U-Atyp','1E8X-Atyp','1IA9-Atyp','1B6C-TKL','1F3M-STE','1CSN-CK1']
#x2=load('/home/linc/c00219805/Research/Protien_Database/S2/jaccard_Gen_CE_localFeatures_scop.txt');
#z2=linkage(x2,'average');
#h=figure;
#dendrogram(z2,31,'label',fileList,'orientation','left');
#saveas(h,'avg_ce-local-dendrogram_NoLeafOrder_Gen','jpg');

input_entity_files = '/home/linc/c00219805/Research/Protien_Database/S2/'

def get_hierarchical_clustering(dist,X):
	linkage_matrix = average(dist)   
	print linkage_matrix
	print linkage_matrix.shape, dist.shape

	#Checking the Cophenetic Correclation Coefficient to compare 
	#actual pairwise distances with HClustering implied distance.
	c, coph_dists = cophenet(linkage_matrix, pdist(X))
	print 'Cophenetic correlation: ', c #This gave 0.336
	knee = np.diff(linkage_matrix[::-1, 2], 2)
	print knee,knee.argmax()
	print hac.fcluster(linkage_matrix,1)
	fig, ax = plt.subplots(figsize=(15, 20)) # set size
	ax = dendrogram(linkage_matrix, orientation="left", labels=fileList) 
	plt.tick_params(\
	    axis= 'x',          # changes apply to the x-axis
	    which='both',      # both major and minor ticks are affected
	    )
	plt.axvline(x=0.85, c='k')
	#plt.tight_layout() #show plot with tight layout
	plt.savefig(input_entity_files +'dendo'+setting+'.png', dpi=200)
	plt.show()
def get_hierarchical_clustering_III(X,xs,ys):
	"""Using the method Joenhees mentioned."""
	# generate the linkage matrix
	Z = linkage(X, 'average')
	c, coph_dists = cophenet(Z, pdist(X))
	print 'Cophenetic correlation: ', c #This gave 0.875
	print Z,X.shape
	print X[[8,17,18]]
	idxs = [8,17,18]
	plt.figure(figsize=(10, 8))
	#plt.scatter(X[:,0], X[:,1])  # plot all points
	#plt.scatter(xs[idxs], ys[idxs], c='r')  # plot interesting points in red again
	#plt.show()
	plt.scatter(xs, ys,label=fileList)
	plt.scatter(xs[idxs], ys[idxs], c='r')
	idxs = [2]
	plt.scatter(xs[idxs], ys[idxs], c='y')
	for i, txt in enumerate(fileList):
		plt.annotate(txt, (xs[i],ys[i]))
	#plt.show()

	plt.figure(figsize=(25, 10))
	plt.title('Hierarchical Clustering Dendrogram')
	plt.xlabel('sample index')
	plt.ylabel('distance')
	dendrogram(
	    Z,
	    leaf_rotation=90.,  # rotates the x axis labels
	    leaf_font_size=8.,  # font size for the x axis labels
	    labels=fileList,
	   
	)
	plt.axhline(y=15800, c='k')
	plt.savefig(input_entity_files +'Hierarchical_clustering_dendo_'+setting+'.png', dpi=200)
	plt.show()

	#clusters=hac.fcluster(Z,15800,criterion='distance')
	#print clusters

	clusters=hac.fcluster(Z,8,criterion='maxclust')
	#print clusters

	plt.figure(figsize=(10, 8))
	plt.scatter(xs, ys, c=clusters, cmap='prism')  # plot points with cluster dependent colors
	for i, txt in enumerate(fileList):
		plt.annotate(txt, (xs[i],ys[i]))
	plt.savefig(input_entity_files +'Hierarchical_clustering'+setting+'.png', dpi=200)
	plt.show()

	#Number of clusters
	last = Z[-30:, 2]
	last_rev = last[::-1]
	idxs = np.arange(1, len(last) + 1)
	plt.plot(idxs, last_rev)

	acceleration = np.diff(last, 2)  # 2nd derivative of the distances
	acceleration_rev = acceleration[::-1]
	plt.plot(idxs[:-2] + 1, acceleration_rev)
	#plt.show()
	k = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters
	print "clusters:", k

def get_hierarchical_clustering_II(cx,no_of_clusters,fileList):
	
	Hclustering = AgglomerativeClustering(n_clusters=no_of_clusters, affinity='manhattan', linkage='average')
	Hclustering.fit(cx)
	print Hclustering.labels_
	ms = np.column_stack((fileList,Hclustering.labels_))
	df = pd.DataFrame(ms,columns = ['Ground truth','Clusters'])
	print pd.crosstab(df['Ground truth'], df['Clusters'],margins=True)

def get_kmeans_clustering(cx,no_of_clusters,fileList):
	num_clusters = no_of_clusters
	km = KMeans(n_clusters=num_clusters)
	km.fit(cx)
	clusters = km.labels_.tolist()
	frame = pd.DataFrame({'protein': fileList, 'cluster':clusters}, index = [clusters] , columns = ['protein', 'cluster'])
	#print frame.sort_values(by = 'cluster')
	return km,frame
def cluster_analysis(frame,no_of_clusters,xs, ys):
	clusters = frame['cluster'].values#km.labels_.tolist()
	colors = ['#1b9e77',  '#d95f02',  '#7570b3', '#e7298a', '#66a61e', '#00FFFF',  '#000080',  '#00FF00',  '#FFFF00',  '#808080','#DAF7A6','#581845','#EBDEF0','#F9EBEA','#EAFAF1']
	
	clus =['Cluster 0','Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Cluster 7','Cluster 8','Cluster 9','Cluster 10','Cluster 11','Cluster 12','Cluster 13','Cluster 14']
	cluster_names = dict(zip(sorted(frame['cluster'].unique()),clus))
	cluster_colors = dict(zip(sorted(frame['cluster'].unique()),colors))
	print cluster_names, cluster_colors
                
	#print frame['cluster'].value_counts()
	#order_centroids = km.cluster_centers_.argsort()[:, ::-1] 
	#print order_centroids
	#for i in range(no_of_clusters):
	#	print order_centroids[i,:6]
	
	df = pd.DataFrame(dict(x=xs, y=ys, label=clusters, title=frame['protein'].values)) 
	print df

	#group by cluster
	groups = df.groupby('label')
	print frame['cluster'].unique()
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

    
	plt.savefig(input_entity_files +'Kmeans_clustering'+setting+setting2+'.png', dpi=200)
	plt.show() #show the plot

def get_two_way_clustering(cx):

	num_clusters = 9
	km = KMeans(n_clusters=num_clusters)
	km.fit(cx)
	clusters = km.labels_.tolist()
	Kx = km.cluster_centers_
	Kx_mapping = {case:cluster for case,cluster in enumerate(km.labels_)}
	print Kx_mapping
	Hclustering = AgglomerativeClustering(n_clusters=9,affinity='cosine', linkage='complete')
	Hclustering.fit(Kx)
	H_mapping = {case:cluster for case,cluster in enumerate(Hclustering.labels_)}
	final_mapping = {case:H_mapping[Kx_mapping[case]] for case in Kx_mapping}
	ms = np.column_stack((fileList, [final_mapping[n] for n in range(max(final_mapping)+1)]))
	df = pd.DataFrame(ms,columns = ['Ground truth','Clusters'])
	print pd.crosstab(df['Ground truth'], df['Clusters'],margins=True)
 
def get_two_way_clustering_II(cx):
	fileList=['1HCK-CMGC','1GNG-CMGC','1HOW-CMGC','1JNK-CMGC','1LUF-TK','2SRC-TK','1M14-TK','1IR3-TK','1GJO-TK','1O6L-AGC','1H1W-AGC','1CDK-AGC','1OMW-AGC','1A06-CAMK','1KWP-CAMK','1TKI-CAMK','1JKL-CAMK','1PHK-CAMK','1IA8-CAMK','1O6Y-Other','1MUO-Other','1LP4-Other','1CJA-Atyp','1BO1-Atyp','1NW1-Atyp','1J7U-Atyp','1E8X-Atyp','1IA9-Atyp','1B6C-TKL','1F3M-STE','1CSN-CK1']
	km,frame = get_kmeans_clustering(cx,9,fileList)
	print Counter(km.labels_.tolist()), Counter(km.labels_.tolist()).most_common(1)[0][0]
	#Performing KMEANs again on the highest points cluster to break it down further
	indeces = np.where(km.labels_ == Counter(km.labels_.tolist()).most_common(1)[0][0])[0]
	#frame = frame.loc[frame['cluster']!=Counter(km.labels_.tolist()).most_common(1)[0][0]]
	print frame
	cx = cx[[indeces],:]
	fileList = [fileList[x] for x in indeces]
	km, frame2 = get_kmeans_clustering(cx[0],7,fileList)
	frame2['cluster'] = frame2['cluster'] +10
	frame_merged = frame.merge(frame2,on=['protein'], how='left')
	print frame_merged
	frame_merged['cluster_y'] = frame_merged['cluster_y'].fillna(frame_merged['cluster_x'])
	frame_merged.apply(pd.to_numeric, errors='ignore')
	frame_merged = frame_merged.drop('cluster_x', 1)
	frame_merged.columns = ['protein', 'cluster']
	print frame_merged
	return km, frame_merged
 
def get_two_way_clustering_III(cx):
	fileList=['1HCK-CMGC','1GNG-CMGC','1HOW-CMGC','1JNK-CMGC','1LUF-TK','2SRC-TK','1M14-TK','1IR3-TK','1GJO-TK','1O6L-AGC','1H1W-AGC','1CDK-AGC','1OMW-AGC','1A06-CAMK','1KWP-CAMK','1TKI-CAMK','1JKL-CAMK','1PHK-CAMK','1IA8-CAMK','1O6Y-Other','1MUO-Other','1LP4-Other','1CJA-Atyp','1BO1-Atyp','1NW1-Atyp','1J7U-Atyp','1E8X-Atyp','1IA9-Atyp','1B6C-TKL','1F3M-STE','1CSN-CK1']
	km,frame = get_kmeans_clustering(cx,9,fileList)
	print Counter(km.labels_.tolist())
	indeces = np.where(km.labels_ == Counter(km.labels_.tolist()).most_common(1)[0][0])[0]
	cx = cx[[indeces],:]
	fileList = [fileList[x] for x in indeces]
	get_hierarchical_clustering_II(cx[0],7,fileList)
	
#Calculating Jaccard similarity as distance measure. X2 is Jaccard Similarity

#x2 = np.loadtxt(open(input_entity_files +"jaccard_Gen_CE_localFeatures_scop_titlinew.txt", "rb"), delimiter=",")
reader = csv.reader(open(input_entity_files +"jaccard_Gen_CE_localFeatures_scop"+setting+".txt", "rb"), delimiter=",")
x = list(reader)
j = 0
for i in x:
	j +=1
x2 = np.array(x).astype("float")

#This is raw matrix, without distance calculations
cx = np.loadtxt(open(input_entity_files +"localFeatureVect"+ setting +".csv", "rb"), delimiter=",")
print cx.shape
doc_clean = [str(i) for i in range(0,149349)]
#print dict(zip(sorted([i for i in range(0,149349)]),doc_clean))
dictionary = corpora.Dictionary([doc_clean])
print dictionary

# Creating the object for LDA model using gensim library
Lda = gensim.models.ldamodel.LdaModel

# Running and Trainign LDA model on the document term matrix.
ldamodel = Lda(cx, num_topics=3, id2word = dict(zip(sorted([i for i in range(0,149349)]),doc_clean)), passes=50)

print ldamodel.print_topics(num_topics=3, num_words=3)

print y


#Calculating Cosine similarity as distance measure
cosine = 1 - cosine_similarity(cx)
#Calculating Euclidean as distance measure
euclidean = euclidean_distances(cx,cx)

#Converts multi dimensional data to two dimensions, uses distance matrix
MDS()
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=1)
#use x2 for jaccard or cosine for cosine
pos = mds.fit_transform(euclidean)
xs, ys = pos[:, 0], pos[:, 1]
print pd.DataFrame(dict(x=xs, y=ys))


#print '---------------------------------------------Method 1-----------------------------------------------------'
#Hierarchical
get_hierarchical_clustering(x2,cx)
#get_hierarchical_clustering_III(cx,xs, ys)



#print '---------------------------------------------Method 2-----------------------------------------------------'
#KMEANS
#Single clustering
#km, frame = get_kmeans_clustering(cx,9,fileList)

#Two level clustering
#km, frame =get_two_way_clustering_II(cx)

#Cluster Analysis
cluster_analysis(frame,9,xs, ys)

#get_two_way_clustering_III(cx)
