import glob, re, pandas as pd

def test():
	path = '/home/linc/c00219805/Research/Protien_Database/'
	subfolder= '/extracted_new_samples/sample_wu//'
	fileType = "*.pdb"       
	files=glob.glob(path+subfolder+fileType) 
	count = 0  
	chains = [] 
	for fileName in files:
		chains = [] 
		inFile=open(fileName,'r')
		for i in inFile:
		    if (i[0:6].rstrip()=="NUMMDL"):
		                numOfModels=i[10:14].rstrip()
		    if ((i[0:6].rstrip()=="ENDMDL") ):
		                break
		    if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
		                break
		               
		    if(i[0:4].rstrip())=="ATOM" and (i[13:15].rstrip())=="CA" and i[17:20]!= "UNK" :
		    	count += 1
		    	chains.append(i[21])
		    	#print count
		    	#print i[16], i[13:15],i[21]
		print fileName,set(chains)
	    	
		    # else:
		    # 	print i
		    # 	print 'False'
	req_lines = []
	for fileName in files:
	    	with open(fileName) as f:
	    	    lines = f.readlines()
	    	    for line in lines:
	    	    	if ('ATOM' in line) and ('REMARK' not in line) and ('UNK' not in line) :
	    	    		line = line.strip()
	    		        # split each column on whitespace
	    		        req_lines.append(re.split('\s+', line, maxsplit=10))
	df = pd.DataFrame(req_lines)
	df.columns = ['type','id','atom','aminoacid','chain','aaid','xcord','ycord','zcord','a','b']
	#df['sortedid']= df.apply(lambda row: atomsCode[row['atom']],axis=1)
	#df['id'] = df[['id']].apply(pd.to_numeric)
	#df['aaid'] = df[['aaid']].apply(pd.to_numeric)
	#print df

def copy_folders():
	a = [25,26,27,37,40]
	bins = range(min(a), max(a)+2,2)
	 
	print pd.cut(a, bins,include_lowest = True)
	
copy_folders()