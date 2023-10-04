import math
import cmath
import sys
import random
import glob
import time
import os,os.path
from scipy import spatial

class SimilarityScore:

    def __init__(self,**kwargs):
        self.setting =  kwargs["setting"]
        self.f2_out = open(kwargs["outFolder"]+'//jaccard_similarity'+self.setting+'.txt','w')
        self.f2_out_cosine = open(kwargs["outFolder"]+'//cosine_similarity'+self.setting+'.txt','w')
        self.f2_out_dice = open(kwargs["outFolder"]+'//dice_similarity'+self.setting+'.txt','w')
        self.fileNames = open(kwargs["outFolder"]+'//localFeatureVect'+self.setting +'.csv','r')
        self.n=len(kwargs["filesList"])
        self.f = kwargs["filesList"]

    def calculate_similarity(self):
        num=0
        dict_={}
        dict_freq = {}
        for i in self.fileNames:
            i=i.split(',')
            dict2={}
            numOfCols=len(i)
            max_norm =  max([float(p) for p in i])
            for k in range(numOfCols-1):#leave the last entry out as it is the dummy 0 entry for all instances
                dict2[k]=float(i[k].rstrip())/max_norm
                #print dict2[k],float(i[k].rstrip()),max_norm
                if k in dict_freq:
                    dict_freq[k].append(float(i[k].rstrip())/max_norm)
                else:
                    dict_freq[k]=[]
                    dict_freq[k].append(float(i[k].rstrip())/max_norm)
            dict_[num]=dict2
            num+=1
        
        for i in range(num):

            col=len(dict_[i])
            for j in range(num):
                print i,j
                result = 1 - spatial.distance.cosine(dict_[i].values(), dict_[j].values())
                result_jaccard = 1 - spatial.distance.jaccard(dict_[i].values(), dict_[j].values())
                result_dice = 1 - spatial.distance.dice(dict_[i].values(), dict_[j].values())
                if j == num-1:
                        #To eliminate last comma
                        self.f2_out.writelines([str(result_jaccard)])
                        self.f2_out_cosine.writelines([str(result)])
                        self.f2_out_dice.writelines([str(result_dice)])
                else:
                        self.f2_out.writelines([str(result_jaccard),','])
                        self.f2_out_cosine.writelines([str(result),','])
                        self.f2_out_dice.writelines([str(result_dice),','])

            self.f2_out.writelines(['\n'])
            self.f2_out_cosine.writelines(['\n'])
            self.f2_out_dice.writelines(['\n'])
        self.f2_out.close()
        self.f2_out_cosine.close()
        self.f2_out_dice.close()
            

