import time
import os,os.path
import pandas as pd
from scipy import spatial
from itertools import combinations
import numpy as np
import StringIO
from joblib import Parallel, delayed, cpu_count

# def unwrap_self(arg, **kwarg):
#     print arg[0]
#     return JaccardCoefficient.process_files(arg[0], arg[1])

class JaccardCoefficient:

    def __init__(self,**kwargs):
        self.setting =  kwargs["setting"]
        
        #self.files = open(kwargs["outFolder"]+'//localFeatureVect'+self.setting +'.csv','r')
        self.files = kwargs["outFolder"]+'//localFeatureVect'+self.setting +'.csv'
        #self.allSimilarityCSV = kwargs["outFolder"]+'//similarity_measures_and_values'+self.setting+'.csv'
        #self.n=len(kwargs["filesList"])
        #self.f = kwargs["filesList"]
        #self.normalize = kwargs["normalize"]
        self.sample_dict = kwargs["sample_dict"]

        self.f2_out_normal = kwargs["outFolder"]+'//normal_jaccard_similarity'+self.setting+'.csv'
        self.f2_out_generalised = kwargs["outFolder"]+'//generalised_jaccard_similarity'+self.setting+'.csv'
        self.f2_out_wu = kwargs["outFolder"]+'//wu_jaccard_similarity'+self.setting+'.csv'
        self.f2_out_sarika1 = kwargs["outFolder"]+'//sarika_jaccard1_similarity'+self.setting+'.csv'
        self.f2_out_cosine = kwargs["outFolder"]+'//cosine_similarity'+self.setting+'.csv'


    def process_files(self):
            arrs = []
            filenames = []
            print self.files
            with open(self.files) as fcsv:
                lines=fcsv.readlines()
                #print lines
                for idx,line in enumerate(lines):
                    filenames.append(self.sample_dict[str(line.split(';')[0]).upper()] \
                        + '-' + str(line.split(';')[0]).upper() )
                    #filenames.append(line.split(';')[0])
                    l = list(line.split(';')[1].split(','))
                    l_arr = np.asarray(l[:-1]).astype(np.float) 
                    arrs.append(l_arr)
            data = np.array(arrs)
            #print data
            data_sum = np.sum(data,axis=1)
            data_jac = np.copy(data)
            data_jac[data_jac>0]=1

            lst_a = np.arange(data.shape[0])

            lst_cmb = list(combinations(lst_a,2))
            #print(lst_cmb)

            normal = np.zeros((data.shape[0],data.shape[0]))
            generalised = np.zeros_like(normal)
            sarika = np.zeros_like(normal)
            wu = np.zeros_like(normal)
            cosine = np.zeros_like(normal)
            #name = i.split(';')[0]
            #print('file:{}'.format(i.split(';')[0]))
            #self.fileNames.append(self.sample_dict[str(i.split(';')[0]).upper()]+ '-' + str(i.split(';')[0]).upper() )
            
            for c in lst_cmb:
                idx_a, idx_b = c
                a = data[idx_a]
                a_sum = data_sum[idx_a]
                a_jac = data_jac[idx_a]
                b = data[idx_b]
                b_sum = data_sum[idx_b]
                b_jac = data_jac[idx_b]

                non_zeros = (a >0) & (b > 0)
                summed_array = a + b

                numerator_jac = np.sum(np.minimum(a_jac,b_jac))
                denomenator_jac = np.sum(np.maximum(a_jac,b_jac))
                numerator_gen_jac =np.sum(np.minimum(a,b))
                denomenator_gen_jac =np.sum(np.maximum(a,b))
                num_sim = np.sum(summed_array[non_zeros])
                #result = 1 - spatial.distance.cosine(a, b)
                result = 1 #- spatial.distance.cosine(a, b)

                if (denomenator_jac == 0):
                    print('There is something wrong. Denominator is Zero! ', idx_a, idx_b, numerator_jac, denomenator_jac)
                else:
                    dist_gen_jac=1.0-(float(numerator_gen_jac)/float(denomenator_gen_jac))                    
                    dist_jac=1.0-(float(numerator_jac)/float(denomenator_jac))

                    denomenator_wu = min(denomenator_gen_jac,max(a_sum,b_sum) )
                    dist_wu = 1.0-(float(numerator_gen_jac)/float(denomenator_wu))
                    
                    numerator_sarika = num_sim
                    denomenator_sarika = a_sum+b_sum
                    dist_sarika = 1.0-(float(numerator_sarika)/float(denomenator_sarika))

                normal[idx_a,idx_b] = dist_jac
                normal[idx_b,idx_a] = dist_jac
                generalised[idx_a,idx_b] = dist_gen_jac
                generalised[idx_b,idx_a] = dist_gen_jac
                sarika[idx_a,idx_b] = dist_sarika
                sarika[idx_b,idx_a] = dist_sarika
                wu[idx_a,idx_b] = dist_wu
                wu[idx_b,idx_a] = dist_wu
                cosine[idx_a,idx_b] = result*100
                cosine[idx_b,idx_a] = result*100
            pd.DataFrame(normal, columns = filenames).to_csv(self.f2_out_normal)
            pd.DataFrame(generalised, columns = filenames).to_csv(self.f2_out_generalised)
            pd.DataFrame(sarika, columns = filenames).to_csv(self.f2_out_sarika1)
            pd.DataFrame(wu, columns = filenames).to_csv(self.f2_out_wu)
            pd.DataFrame(cosine, columns = filenames).to_csv(self.f2_out_cosine)

    def calculate_jaccard(self):
        
        start_time=time.time()
        
        # results   = Parallel(n_jobs=cpu_count() - 1, verbose=10, \
        #                 backend="multiprocessing", batch_size="auto")(delayed(unwrap_self)(i) for i in zip([self]*len(lines),lines))

        self.process_files()
            
        end_time=time.time()
        total_time=((end_time)-(start_time))
        print("Time taken for writing to files: {}".format(total_time))

        
        

