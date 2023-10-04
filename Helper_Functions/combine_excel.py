import os
import pandas as pd
import xlrd
import argparse

__author__ = "Venkata Sarika Kondra"
__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"

parser = argparse.ArgumentParser(description='Search Keys.')
parser.add_argument('--input', '-input', metavar='input', \
	default='C:\Users\wxx6941\Desktop\freq_less5_dist_less6_triplets.xlsx', \
	help='Name of the input file.')
parser.add_argument('--output', '-output', metavar='output', \
	default='C:\\Users\\wxx6941\\Desktop\\single_freq_less5_dist_less6_triplets.csv', \
	help='Name of the output file.')
	

args = parser.parse_args()
column_names = [u'Unnamed: 0', u'key', u'aa0', u'pos0', u'aa1', u'pos1', u'aa2',\
       u'pos2', u'classT', u'theta', u'classL', u'distance', u'x0', u'y0',\
       u'z0', u'x1', u'y1', u'z1', u'x2', u'y2', u'z2']
sheets = pd.ExcelFile(args.input).sheet_names
df_all = pd.DataFrame(columns = column_names)
xl = pd.ExcelFile(args.input)
for sheet in xl.sheet_names:
	df = xl.parse(sheet, column = column_names)
	df['fileName'] = sheet
	df_all = df_all.append(df)
print df_all
df_all.to_csv(args.output)