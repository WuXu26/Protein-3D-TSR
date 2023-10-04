import copy
import os,os.path

class Vectorization:

    def __init__(self, **kwargs):
        self.setting = kwargs["setting"] 
        self.fileList = kwargs["filesList"]
        self.f1_in = open(kwargs["outFolder"] + '//localFeatureSelection' + self.setting + '.txt', 'r')
        self.f1_out = open(kwargs["outFolder"] + '//localFeatureVect' + self.setting + '.csv', 'w')

    def vectorize(self):
        keyDict = {}
        for i in self.f1_in:
            i = i.split()
            keyDict[i[0].rstrip()] = 0
        self.f1_in.close()
        for i in self.fileList:
            f2_in = open(i,'r')
            keyDict1 = copy.deepcopy(keyDict)
            for j in f2_in:
            	j = j.split()
            	if j[0].rstrip() in keyDict1:
            	    keyDict1[j[0].rstrip()] = j[1].rstrip()
            f2_in.close()
            self.f1_out.writelines([str(i.split('.')[0]), ';'])
            for k in keyDict1:
                self.f1_out.writelines([str(keyDict1[k]), ','])
            self.f1_out.writelines(['0', '\n'])
            
        self.f1_out.close()
        
            	
