import numpy as np
import argparse, os, ntpath, glob
from os.path import expanduser
import pandas as pd


parser = argparse.ArgumentParser(description='Parallel Key Generation.')
parser.add_argument('--file_path', '-file_path', metavar='file_path', \
    default=os.path.join(expanduser('~'),'Research', 'Protien_Database', \
    'extracted_new_samples', 'testing','t1','1V2W.pdb'), \
    help='Directory of input sample and other files.')
parser.add_argument('--aacds', '-aacds', metavar='aacds', \
    default='gly,thr',\
    help='Amino acids that are to be searched.')

def get_distance(a,b):
    dist = np.linalg.norm(a-b)
    return dist

if __name__ == '__main__':	
    args = parser.parse_args()
    aa1 = args.aacds.split(',')[0].upper()
    aa2 = args.aacds.split(',')[1].upper()
    df = pd.read_csv(os.path.join(args.file_path, 'sample_details.csv'))
    df_dict = dict(zip(df.protein,df.chain))
    summary = []
    directoryPath = os.path.join(args.file_path)
    for file in os.listdir(directoryPath):
        if(file.endswith(".pdb")):
            print(file[0:4])
            print()
            inFile = open(directoryPath + file,'r')
            result = []
            a_coords = []
            b_coords = []
            #print(inFile, result)
            for i in inFile:
                if (i[0:6].rstrip()=="NUMMDL"):
                    numOfModels=i[10:14].rstrip()
                if ((i[0:6].rstrip()=="ENDMDL")): # or (i[0:6].rstrip()=='TER')):
                    break
                if (i[0:6].rstrip()=="MODEL" and int(i[10:14].rstrip())>1):
                    break
                if(i[0:4].rstrip() == "ATOM"):
                    if(i[16:20].strip() == aa1.upper()):
                        a_coords.append((i[16:20].rstrip(), str(i[22:27]), float(i[30:38]), float(i[38:46]), float(i[46:54])))
                    if(i[16:20].strip() == aa2.upper()) :
                        b_coords.append((i[16:20].rstrip(), str(i[22:27]), float(i[30:38]), float(i[38:46]), float(i[46:54])))
            if aa1 == aa2:
                b_coords = a_coords
            for a in a_coords:
                for b in b_coords:
                    if a[1] != b[1]:
                        a_coord = np.array((a[2], a[3], a[4]))
                        b_coord = np.array((b[2], b[3], b[4]))
                        result.append((a[0], a[1], b[0], b[1], a[2], a[3], a[4], b[2], b[3], b[4], get_distance(a_coord,b_coord)))
             # pd.DataFrame(result, columns = ['aa1', 'aa1_pos','aa2','aa2_pos', 'x_aa1', 'y_aa1', 'z_aa1',  'x_aa2', 'y_aa2', 'z_aa2', 'distance'])\
            # 	.to_excel(writer,sheet_name=ntpath.basename(file)[:4])
            df = pd.DataFrame(result, columns = ['aa1', 'aa1_pos','aa2','aa2_pos', 'x_aa1', 'y_aa1', 'z_aa1',  'x_aa2', 'y_aa2', 'z_aa2', 'distance'])
            #summary.append((ntpath.basename(file)[:4], pd.to_numeric(df['distance']).mean()))
            print(df)
            df.to_csv(file[0:4] + "_Distances_Output")
            #df.to_csv(os.path.join(os.path.dirname(args.file_path),'{}_{}_{}_distances.csv'.format(ntpath.basename(file)[:4], aa1, aa2)))
            inFile.close()
        #pd.DataFrame(summary, columns = ['protein', 'dist_mean']).to_csv(os.path.join(os.path.dirname(args.file_path),'summary_distances.csv'))
            #writer.close()
			
