
import math,numpy,ntpath,glob
import copy,pandas as pd,numpy as np
from pathlib import Path


class FeatureSelection:

    def __init__(self,**kwargs):
      self.setting = kwargs["setting"] 
      self.outfolder = kwargs["outFolder"]
      self.filesList = kwargs["filesList"]
      self.numOfGap = kwargs["numOfGap"]
      self.mad = kwargs["mad"]
      self.keyCombine = kwargs["keyCombine"]
      self.allKeyCountDict={} 
      self.masterDist = {}
      self.allKeysOnePlusDict={} 
      self.generate_dict()
      self.f2_out=open(kwargs["outFolder"] + '//localFeatureSelection' + self.setting + '.txt','w')
      self.fs = kwargs["featureSelection"]

    def generate_dict(self):
        for i in self.filesList:
            localdist = {}
            f1=open(i,'r')
            for j in f1:
                j=j.split()
                key_=j[0].rstrip()
                val_=int(j[1].rstrip())
                localdist[int(key_)] = val_
                if key_ in self.allKeyCountDict:
                    self.allKeyCountDict[key_].append(val_)
                else:
                    self.allKeyCountDict[key_]=[]
                    self.allKeyCountDict[key_].append(val_)
            self.masterDist[i] = localdist
            f1.close()
    def f(self,x):
        return list(filter(lambda a: float(a) != 0, x.split(',')))
    def get_plus_one(self):    
        df_ = pd.DataFrame( columns=['grp','freqList'])
        
        minBin = min(list(map(int,self.allKeyCountDict.keys())))
        maxBin = max(list(map(int,self.allKeyCountDict.keys())))
        bins = range(minBin, maxBin + 2,self.keyCombine + self.keyCombine)   
        labels = range(1, len(bins))
        df_all = pd.DataFrame()
        for file,value in self.masterDist.items():
            if not Path("{}//{}.keys_keycombine{}".format(self.outfolder,ntpath.basename(file)[:4],str(self.keyCombine))).exists():
                df = pd.DataFrame(value.items(),columns = ["key","freq"]).sort_values(by="key")
                df['grp'] = pd.cut(df.key, bins,labels = labels,include_lowest = True)
                new = df[['key','grp']].drop_duplicates().copy()
                df_all = pd.concat([df_all,new])
                #If you file is already not saved to perform pre-processing of common_keys.py
                if not Path("{}//{}.keys_original_keycombine{}".format(self.outfolder,ntpath.basename(file)[:4],str(self.keyCombine))).exists():
                    df.to_csv(self.outfolder + "//" + ntpath.basename(file)[:4] + ".keys_original_keycombine" + str(self.keyCombine), sep='\t',header = True,index=False)
                df = df.drop(["key"],axis = 1)
                grouped = df.groupby([ 'grp'], as_index=False)
                final = grouped.agg(np.sum).dropna()
                final.to_csv("{}//{}.keys_keycombine{}".format(self.outfolder,ntpath.basename(file)[:4],str(self.keyCombine)), sep='\t',header = False,index=False)
        
        changedFiles = glob.glob(self.outfolder + '//*.keys_keycombine' + str(self.keyCombine))
       
        for i in changedFiles:
            f1=open(i,'r')
            for j in f1:
                j=j.split()
                key_=j[0].rstrip()
                val_=float(j[1].rstrip())
                if key_ in self.allKeysOnePlusDict:
                    self.allKeysOnePlusDict[key_].append(val_)
                else:
                    self.allKeysOnePlusDict[key_]=[]
                    self.allKeysOnePlusDict[key_].append(val_)
            f1.close()
         
         
    def feature_selection(self):
        loopThruDict = self.allKeyCountDict
        if  self.keyCombine > 0:
            self.get_plus_one()
            loopThruDict = self.allKeysOnePlusDict 
        
        for keys_ in loopThruDict:
            list_=[]
            list_=map(float,loopThruDict[keys_])
            numOfMatch=len(list_) #number of proteins that have the key
            numOfGap=len(self.filesList)-numOfMatch #number of proteins that DO NOT have the key
            mean_=numpy.mean(list_) #average number of times a key occurs in the list of proteins
            mad_=0.0
            for i in range(len(list_)):
                mad_=mad_ + math.fabs(list_[i]-mean_)
            mad_=float(mad_)/float(numOfMatch)# calculate MAD
            if not self.fs:
                self.f2_out.writelines([str(keys_),'\t',str(mad_),'\n'])
            else:
                if ((numOfGap>=self.numOfGap and mad_ <= self.mad)):# and (numOfMatch>2)):
                    self.f2_out.writelines([str(keys_),'\t',str(mad_),'\n'])
        self.f2_out.close()#Never miss closing files, if you are running a program after this using the file as input it fails.
      
