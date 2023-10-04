import math
import cmath
import sys
import random
import glob
import time
import os,os.path
import pandas as pd
from scipy import spatial
import numpy as np
from numpy import array
from joblib import Parallel, delayed, cpu_count

def unwrap_self(arg, **kwarg):
    print arg[0]
    return JaccardCoefficient.process_files(arg[0], arg[1])

class JaccardCoefficient:

    def __init__(self,**kwargs):
        self.setting =  kwargs["setting"]
        
        self.files = open(kwargs["outFolder"]+'//localFeatureVect'+self.setting +'.csv','r')
        self.allSimilarityCSV = kwargs["outFolder"]+'//similarity_measures_and_values'+self.setting+'.csv'
        self.n=len(kwargs["filesList"])
        self.f = kwargs["filesList"]
        self.normalize = kwargs["normalize"]
        self.sample_dict = kwargs["sample_dict"]

        self.f2_out_normal = kwargs["outFolder"]+'//normal_jaccard_similarity'+self.setting+'.csv'
        self.f2_out_generalised = kwargs["outFolder"]+'//generalised_jaccard_similarity'+self.setting+'.csv'
        self.f2_out_wu = kwargs["outFolder"]+'//wu_jaccard_similarity'+self.setting+'.csv'
        self.f2_out_sarika1 = kwargs["outFolder"]+'//sarika_jaccard1_similarity'+self.setting+'.csv'
        self.f2_out_cosine = kwargs["outFolder"]+'//cosine_similarity'+self.setting+'.csv'

        self.lines = self.files.readlines()
        self.lines = array(self.lines)
        #print type(self.lines), type(array(self.lines))
        self.files.close()
        self.fileNames = []

    def process_files(self,i):
            normal = []
            generalised = []
            sarika = []
            wu = []
            cosine = []
            name = i.split(';')[0]
            print('file:{}'.format(i.split(';')[0]))
            self.fileNames.append(self.sample_dict[str(i.split(';')[0]).upper()]+ '-' + str(i.split(';')[0]).upper() )
            
            i=list(i.split(';')[1].split(','))
            a = np.asarray(i[:-1]).astype(np.float)  
            a_sum = np.sum(a)
            a_jac = np.copy(a)
            a_jac[a_jac>0] = 1
            
            for j in self.lines:
                j = list(j.split(';')[1].split(','))
                b = np.asarray(j[:-1]).astype(np.float)
                non_zeros = (a >0) & (b > 0)
                summed_array = a + b
                b_sum = np.sum(b)
                b_jac = np.copy(b)
                b_jac[b_jac>0] = 1
                
                numerator_jac = np.sum(np.minimum(a_jac,b_jac))
                denomenator_jac = np.sum(np.maximum(a_jac,b_jac))
                numerator_gen_jac =np.sum(np.minimum(a,b))
                denomenator_gen_jac =np.sum(np.maximum(a,b))
                num_sim = np.sum(summed_array[non_zeros])
                result = 1 - spatial.distance.cosine(a, b)
                

                if denomenator_jac == 0:
                        print('There is something wrong. Denominator is Zero! ',i, j, numerator_jac, denomenator_jac)
                else:
                    dist_gen_jac=1.0-(float(numerator_gen_jac)/float(denomenator_gen_jac))                    
                    dist_jac=1.0-(float(numerator_jac)/float(denomenator_jac))

                    denomenator_wu = min(denomenator_gen_jac,max(a_sum,b_sum) )
                    dist_wu = 1.0-(float(numerator_gen_jac)/float(denomenator_wu))
                    
                    numerator_sarika = num_sim
                    denomenator_sarika = a_sum+b_sum
                    dist_sarika = 1.0-(float(numerator_sarika)/float(denomenator_sarika))

                    
                    normal.append(str(dist_jac)) 
                    generalised.append(str(dist_gen_jac))                      
                    sarika.append(str(dist_sarika))
                    wu.append(str(dist_wu))
                    cosine.append(str(result*100))
            return (name, normal, generalised, sarika, wu, cosine)

    def calculate_jaccard(self):
        
        #lines = self.fileNames.readlines()
        #self.fileNames.close()
        lines = self.lines      
        start_time=time.time()
        
        #Add code for normalization later
        # file, normal_all, generalised_all, sarika_all, wu_all, cosine_all   = Parallel(n_jobs=cpu_count() - 1, verbose=10, \
        #                 backend="multiprocessing", batch_size="auto")(\
        #                 delayed(process_files)(i, lines) for i in lines)

        # file, normal_all, generalised_all, sarika_all, wu_all, cosine_all   = Parallel(n_jobs=cpu_count() - 1, verbose=10, \
        #                 backend="multiprocessing", batch_size="auto")(\
        #                 delayed(unwrap_self)(i, lines) for i in lines)
        #file, normal_all, generalised_all, sarika_all, wu_all, cosine_all
        results   = Parallel(n_jobs=cpu_count() - 1, verbose=10, \
                        backend="multiprocessing", batch_size="auto")(delayed(unwrap_self)(i) for i in zip([self]*len(lines),lines))


        file = [self.sample_dict[str(x[0]).upper()]+ '-' + str(x[0]).upper() for x in results]
        #self.fileNames.append( )
        normal_all = [x[1] for x in results]
        generalised_all = [x[2] for x in results]
        sarika_all = [x[3] for x in results]
        wu_all = [x[4] for x in results]
        cosine_all = [x[5] for x in results]
            
        end_time=time.time()
        total_time=((end_time)-(start_time))
        print("Time taken for writing to files: {}".format(total_time))

        pd.DataFrame(normal_all,columns = file).to_csv(self.f2_out_normal)
        pd.DataFrame(generalised_all,columns = file).to_csv(self.f2_out_generalised)
        pd.DataFrame(sarika_all,columns = file).to_csv(self.f2_out_sarika1)
        pd.DataFrame(wu_all,columns = file).to_csv(self.f2_out_wu)
        pd.DataFrame(cosine_all,columns = file).to_csv(self.f2_out_cosine)
        
        

