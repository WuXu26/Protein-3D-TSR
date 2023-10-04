'''

Code for calculating PDB to TSRs with changed length bins (domain knowledge incorporated)
'''
import os,math,ntpath
import glob
import time

class KeyGeneration:
    #Changed according to Sarika's binning
    def thetaClass_(self, binBoundaries, Theta, type):
        for i in binBoundaries:
            if Theta < binBoundaries[0]:
                if type == 0:
                    print 'out of index'
                else:
                    classL = binBoundaries.index(binBoundaries[0])+1
                break
            if (Theta < i) :
                if type ==0:
                     classL = binBoundaries.index(i) 
                else:
                     classL =binBoundaries.index(i) + 1
                break
        if Theta >= binBoundaries[-1]:
            if type ==0:
                if Theta == binBoundaries[-1]:
                    classL = binBoundaries.index(binBoundaries[-1])
            else :
                 classL = binBoundaries.index(binBoundaries[-1]) +2
        return classL


    def calcDist(self,indexLabel1,indexLabel2):
        x1=xCord[indexLabel1]
        x2=xCord[indexLabel2]
        y1=yCord[indexLabel1]
        y2=yCord[indexLabel2]
        z1=zCord[indexLabel1]
        z2=zCord[indexLabel2]
        distance=(((x1-x2)**2+(y2-y1)**2+(z2-z1)**2)**0.5)
        return distance

    def indexFind(self,index_of_2,i1,j1,k1):
        if index_of_2==i1:
            indexOf0=j1
            indexOf1=k1
        elif index_of_2==j1:
            indexOf0=i1
            indexOf1=k1
        elif index_of_2==k1:
            indexOf0=i1
            indexOf1=j1

        return indexOf0, indexOf1

    def processFiles(self,fileName):
        print fileName
        print self.outFolder +"//"+ntpath.basename(fileName)[:-3]
        start_time=time.time()
        filesDict={}
        inFile=open(fileName,'r')
        outFile2 = open((self.outFolder+"//" +ntpath.basename(fileName)[:-3] + "keys_theta"+str(self.dTheta)+"_dist"+str(self.dLen)), "w")
        print outFile2
        fileTriplets = open((self.outFolder +ntpath.basename(fileName)[:-3] + "triplets_theta"+str(self.dTheta)+"_dist"+str(self.dLen)), "w")

        global xCord, yCord, zCord
        aminoAcidName={}
        xCord={}
        yCord={}
        zCord={}
        seq_number={}
        counter=0

        for i in inFile:
            if (i[0:6].rstrip()=="NUMMDL"):
                numOfModels=i[10:14].rstrip()
            if ((i[0:6].rstrip()=="ENDMDL") or (i[0:6].rstrip()=='TER')):
                break
            if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
                break
                
            if(i[0:4].rstrip())=="ATOM" and (i[13:15].rstrip())=="CA" and (i[16]=='A'or i[16]==' ')and i[17:20]!= "UNK" :
                
                aminoAcidName[counter]=int(self.aminoAcidLabel[i[17:20]])
                xCord[counter]=(float(i[30:38]))
                yCord[counter]=(float(i[38:46]))
                zCord[counter]=(float(i[46:54]))
                seq_number[counter]=str(i[22:27])
                counter+=1
                
        protLen=len(yCord)
        initialLabel=[]
        sortedLabel=[]
        sortedIndex=[]
        outDist={}
        for m in range(0,3):
            initialLabel.append(0)
            sortedLabel.append(0)
            sortedIndex.append(0)
        for i in range(0,protLen-2):
            for j in range(i+1,protLen-1):
                for k in range(j+1, protLen):
                    global i1,j1,k1
                    i1=i
                    j1=j
                    k1=k
                    keepLabelIndex={}
                    keepLabelIndex[aminoAcidName[i]]=i
                    keepLabelIndex[aminoAcidName[j]]=j
                    keepLabelIndex[aminoAcidName[k]]=k
                    initialLabel[0]=aminoAcidName[i]
                    initialLabel[1]=aminoAcidName[j]
                    initialLabel[2]=aminoAcidName[k]
                    sortedLabel=list(initialLabel)
                    sortedLabel.sort(reverse=True)
                    if (sortedLabel[0]==sortedLabel[1]) and (sortedLabel[1]==sortedLabel[2]):
                        dist1_2Temp=self.calcDist(i,j)
                        dist1_3Temp=self.calcDist(i,k)
                        dist2_3Temp=self.calcDist(j,k)
                        if dist1_2Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                            indexOf0=i
                            indexOf1=j
                            indexOf2=k
                        elif dist1_3Temp>=(max(dist1_2Temp,dist1_3Temp,dist2_3Temp)):
                            indexOf0=i
                            indexOf1=k
                            indexOf2=j
                        else:
                            indexOf0=j
                            indexOf1=k
                            indexOf2=i
                    elif(aminoAcidName[i]!=aminoAcidName[j]) and (aminoAcidName[i]!=aminoAcidName[k]) and (aminoAcidName[j]!=aminoAcidName[k]):
                        for index_ in range(0,3):
                            sortedIndex[index_]=keepLabelIndex[sortedLabel[index_]]
                        indexOf0=sortedIndex[0]
                        indexOf1=sortedIndex[1]
                        indexOf2=sortedIndex[2]
                    elif(sortedLabel[0]==sortedLabel[1]) and (sortedLabel[1]!=sortedLabel[2]):
                        indexOf2=keepLabelIndex[sortedLabel[2]]
                        indices=self.indexFind(indexOf2,i,j,k)
                        a=indexOf2
                        b=indices[0]
                        c=indices[1]
                        dist1_3Temp=self.calcDist(b,a)
                        dist2_3Temp=self.calcDist(c,a)
                        if dist1_3Temp>=dist2_3Temp:
                            indexOf0=indices[0]
                            indexOf1=indices[1] 
                        else:
                            indexOf0=indices[1]
                            indexOf1=indices[0]
                    elif(sortedLabel[0]!=sortedLabel[1]) and (sortedLabel[1]==sortedLabel[2]):
                        indexOf0=keepLabelIndex[sortedLabel[0]]
                        indices=self.indexFind(indexOf0,i,j,k)
                        if self.calcDist(indexOf0,indices[0])>= self.calcDist(indexOf0,indices[1]):
                            indexOf1=indices[0]
                            indexOf2=indices[1] 
                        else:
                            indexOf2=indices[0]
                            indexOf1=indices[1]
                    dist01=self.calcDist(indexOf0,indexOf1)
                    s2=dist01/2
                    dist02=self.calcDist(indexOf0,indexOf2)
                    s1=dist02
                    dist12=dist01
                    dist03=self.calcDist(indexOf1,indexOf2)
                    maxDist=max(dist01,dist02,dist03)
                    s3=(((xCord[indexOf0]+xCord[indexOf1])/2-xCord[indexOf2])**2+((yCord[indexOf0]+yCord[indexOf1])/2-yCord[indexOf2])**2+((zCord[indexOf0]+zCord[indexOf1])/2-zCord[indexOf2])**2)**0.5
                    Theta1=180*(math.acos((s1**2-s2**2-s3**2)/(2*s2*s3)))/3.14
                    if Theta1<=90:
                        Theta=Theta1
                    else:
                        Theta=abs(180-Theta1)
                    classT1=self.thetaClass_(self.thetaBounds,Theta,0)
                    classL1=self.thetaClass_(self.distBounds, maxDist,1)

                    ##getting the positions of AminoAcids in sequence
                    position0 = str(seq_number.values()[indexOf0])
                    position1 = str(seq_number.values()[indexOf1])
                    position2 = str(seq_number.values()[indexOf2])

                    aacd0 = self.aminoAcidLabel.keys()[self.aminoAcidLabel.values().index(aminoAcidName[indexOf0])]
                    aacd1 = self.aminoAcidLabel.keys()[self.aminoAcidLabel.values().index(aminoAcidName[indexOf1])]
                    aacd2 = self.aminoAcidLabel.keys()[self.aminoAcidLabel.values().index(aminoAcidName[indexOf2])]

                    x0 = str(xCord.get(indexOf0))
                    y0 = str(yCord.get(indexOf0))
                    z0 = str(zCord.get(indexOf0))

                    x1 = str(xCord.get(indexOf1))
                    y1 = str(yCord.get(indexOf1))
                    z1 = str(zCord.get(indexOf1))

                    x2 = str(xCord.get(indexOf2))
                    y2 = str(yCord.get(indexOf2))
                    z2 = str(zCord.get(indexOf2))

                    key_2=self.dLen*self.dTheta*(self.numOfLabels**2)*(aminoAcidName[indexOf0]-1)+self.dLen*self.dTheta*(self.numOfLabels)*(aminoAcidName[indexOf1]-1)+self.dLen*self.dTheta*(aminoAcidName[indexOf2]-1)+self.dTheta*(classL1-1)+(classT1-1)
                   
                    if key_2 in filesDict:
                        filesDict[key_2]+=1
                    else:
                        filesDict[key_2]=1

                    line = (str(key_2)+"\t"+str(aacd0)+"\t"+str(position0)+"\t"+str(aacd1)+"\t"+str(position1)+"\t"+str(aacd2)+"\t"+str(position2)+"\t"+str(classT1)+"\t"+str(Theta)+"\t"+str(classL1)+"\t"+str(maxDist)+"\t"+x0+"\t"+y0+"\t"+z0+"\t"+x1+"\t"+y1+"\t"+z1+"\t"+x2+"\t"+y2+"\t"+z2+"\n")
                    fileTriplets.writelines(line)

        for value_ in filesDict:
            outFile2.writelines([str(value_),'\t', str(filesDict[value_]),'\n'])

        end_time=time.time()
        total_time=((end_time)-(start_time))
        print ("FILENAME=",fileName,'\t',"TIME IN SEC(Start,End,Total)=",start_time,',',end_time,',',total_time,'\t',"NUM OF AMINOACIDS=",protLen)

    def __init__(self,**kwargs):
        self.thetaBounds = kwargs["thetaBounds"]
        self.dTheta = len(self.thetaBounds) - 1
        self.distBounds =  kwargs["distBounds"]
        self.dLen = len(self.distBounds) + 1
        self.numOfLabels = kwargs["numOfLabels"]
        self.path = kwargs["path"]
        self.subfolder= kwargs["subFolder"]
        self.aminoAcidCode=kwargs["aminoAcidCode"]
        self.aminoAcidLabel={}
        #aminoAcidGroup={}
        for amino in self.aminoAcidCode:
            amino=amino.split()
            self.aminoAcidLabel[amino[0]]=int(amino[1])
        self.aminoAcidCode.close()
        self.outFolder = self.path+self.subfolder+"theta"+str(self.dTheta)+"_dist"+str(self.dLen)
        #Create output directory
        if not os.path.exists(self.outFolder):
            os.makedirs(self.outFolder)

        
