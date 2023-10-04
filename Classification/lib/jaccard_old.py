import math
import cmath
import sys
import random
import glob
import time
import os,os.path
import pandas as pd
from scipy import spatial

def unwrap_self(arg, **kwarg):
    return square_class.square_int(*arg, **kwarg)

class JaccardCoefficient:

    def __init__(self,**kwargs):
        self.setting =  kwargs["setting"]
        self.f2_out_normal = open(kwargs["outFolder"]+'//normal_jaccard_similarity'+self.setting+'.txt','w')
        self.f2_out_generalised = open(kwargs["outFolder"]+'//generalised_jaccard_similarity'+self.setting+'.txt','w')
        self.f2_out_wu = open(kwargs["outFolder"]+'//wu_jaccard_similarity'+self.setting+'.txt','w')
        self.f2_out_sarika1 = open(kwargs["outFolder"]+'//sarika_jaccard1_similarity'+self.setting+'.txt','w')
        #self.f2_out_sarika2 = open(kwargs["outFolder"]+'//sarika_jaccard2_similarity'+self.setting+'.txt','w')
        self.f2_out_cosine = open(kwargs["outFolder"]+'//cosine_similarity'+self.setting+'.txt','w')
        self.fileNames = open(kwargs["outFolder"]+'//localFeatureVect'+self.setting +'.csv','r')
        self.allSimilarityCSV = kwargs["outFolder"]+'//similarity_measures_and_values'+self.setting+'.csv'
        self.n=len(kwargs["filesList"])
        self.f = kwargs["filesList"]
        self.normalize = kwargs["normalize"]

    

    def calculate_jaccard(self):
        
        num=0
        dict_={}
        dict_freq = {}
        fileNames = []
        start_time=time.time()
        for i in self.fileNames:
            fileNames.append(i.split(';')[0])
            i=i.split(';')[1].split(',')
            dict2={}
            numOfCols=len(i)
            if self.normalize:
                max_norm =  max([float(p) for p in i])
            for k in range(numOfCols-1):#leave the last entry out as it is the dummy 0 entry for all instances
                dict2[k]=float(i[k].rstrip())

                if self.normalize:
                    dict2[k]=float(i[k].rstrip())/max_norm
                    if k in dict_freq:
                        dict_freq[k].append(float(i[k].rstrip())/max_norm)
                    else:
                        dict_freq[k]=[]
                        dict_freq[k].append(float(i[k].rstrip())/max_norm)
            dict_[num]=dict2
            num+=1
        self.fileNames.close()

        end_time=time.time()
        total_time=((end_time)-(start_time))
        print("Time taken for dictionary creation: {}".format(total_time))
        sim_measure = []

        for i in range(num):
            print(i,self.f[i].split('.')[0],fileNames[i])
            self.f2_out_normal.writelines([str(fileNames[i]),','])
            col=len(dict_[i])
            for j in range(num):
                #print self.f[j].split('.')[0]
                result = 1 - spatial.distance.cosine(dict_[i].values(), dict_[j].values())
                numerator_jac=0
                denomenator_jac=0
                numerator_gen_jac=0
                denomenator_gen_jac=0
                a_sum =0
                b_sum =0
                num_sim = 0
                for k in range(col):
                    
                    a=dict_[i][k]
                    b=dict_[j][k] 
                    
                    if self.normalize:
                        if a != 0:
                            a=max(dict_freq[k]) - dict_[i][k]
                        if b != 0:
                            b=max(dict_freq[k]) -dict_[j][k]
                    if a>0:
                        a_jac=1
                    else:
                        a_jac=0
                    if b>0:
                        b_jac=1
                    else:
                        b_jac=0
                    numerator_jac+=min(a_jac,b_jac)
                    denomenator_jac+=max(a_jac,b_jac)
                    numerator_gen_jac+=min(a,b)
                    denomenator_gen_jac+=max(a,b)
                    a_sum+=a
                    b_sum+=b
                    if min(a,b) != 0:
                        num_sim += a+b

                if denomenator_jac==0:
                        print('There is something wrong. Denominator is Zero! ',i, j, numerator_jac, denomenator_jac)
                else:
                    #print('Generalized Jaccard ----numerator:',float(numerator_gen_jac), 'denominator:',float(denomenator_gen_jac))
                    dist_gen_jac=1.0-(float(numerator_gen_jac)/float(denomenator_gen_jac))
                    #dist_gen_jac=100*(float(numerator_gen_jac)/float(denomenator_gen_jac))
                    #sim_measure.append( ['Generalized Jaccard',self.f[i].split('.')[0], self.f[j].split('.')[0],
                     #   numerator_gen_jac,denomenator_gen_jac,1.0-dist_gen_jac,dist_gen_jac])

                    #print('Normal----numerator', 'denominator',float(numerator_jac),float(denomenator_jac))
                    dist_jac=1.0-(float(numerator_jac)/float(denomenator_jac))
                    #dist_jac=100*(float(numerator_jac)/float(denomenator_jac))
                    #sim_measure.append( ['Normal Jaccard',self.f[i].split('.')[0], self.f[j].split('.')[0],
                     #   numerator_jac,denomenator_jac,1.0-dist_jac,dist_jac])

                    denomenator_wu = min(denomenator_gen_jac,max(a_sum,b_sum) )
                    #print('Modified Generalized(Wu)----numerator', 'denominator',float(numerator_gen_jac),float(denomenator_wu))
                    dist_wu = 1.0-(float(numerator_gen_jac)/float(denomenator_wu))
                    #dist_wu = 100*(float(numerator_gen_jac)/float(denomenator_wu))
                    #sim_measure.append( ['Modified Generalized(Wu)',self.f[i].split('.')[0], self.f[j].split('.')[0],
                     #   numerator_gen_jac,denomenator_wu,1.0-dist_wu,dist_wu])

                    numerator_sarika = num_sim
                    denomenator_sarika = a_sum+b_sum
                    #print('Modified Generalized(Sarika)----numerator', 'denominator',float(numerator_sarika),float(denomenator_sarika))
                    dist_sarika = 1.0-(float(numerator_sarika)/float(denomenator_sarika))
                    #dist_sarika = 100*(float(numerator_sarika)/float(denomenator_sarika))
                    #sim_measure.append( ['Modified Generalized(Sarika)',self.f[i].split('.')[0], self.f[j].split('.')[0],
                     #   numerator_sarika,denomenator_sarika,1.0 - dist_sarika,dist_sarika])
                    
                    #To eliminate last comma
                    if j == num-1:                        
                        self.f2_out_normal.writelines([str(dist_jac)])
                        self.f2_out_generalised.writelines([str(dist_gen_jac)])
                        self.f2_out_wu.writelines([str(dist_jac)])
                        self.f2_out_sarika1.writelines([str(dist_sarika)])
                        self.f2_out_cosine.writelines([str(result*100)])
                    else:
                        self.f2_out_normal.writelines([str(dist_jac),','])
                        self.f2_out_generalised.writelines([str(dist_gen_jac),','])
                        self.f2_out_wu.writelines([str(dist_jac),','])
                        self.f2_out_sarika1.writelines([str(dist_sarika),','])
                        self.f2_out_cosine.writelines([str(result*100),','])
            self.f2_out_normal.writelines(['\n'])
            self.f2_out_generalised.writelines(['\n'])
            self.f2_out_wu.writelines(['\n'])
            self.f2_out_sarika1.writelines(['\n'])
            self.f2_out_cosine.writelines(['\n'])
        #result_df = pd.DataFrame(sim_measure,\
         #       columns=['Similarity Measure Used','Protein1','Protein2','numerator','denominator','similarity (num/denom)','dissimilarity(1-(num/denom))'])
        #result_df = result_df[result_df['Protein1'] != result_df['Protein2']]
        #print(result_df)
        #result_df.to_csv(self.allSimilarityCSV)
        self.f2_out_normal.close()
        self.f2_out_generalised.close()
        self.f2_out_wu.close()
        self.f2_out_sarika1.close()
        self.f2_out_cosine.close()
            

