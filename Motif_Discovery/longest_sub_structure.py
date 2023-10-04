#!/usr/bin/env python
#longest_sub_structure.py
"""
	Searches for specific keys and retrieves their details from triplets files
"""

import glob, os, csv, ntpath,socket,argparse, time
import pandas as pd, numpy as np
from collections import Counter


__author__ = "Venkata Sarika Kondra"

__version__ = "1.0.1"
__maintainer__ = "Venkata Sarika Kondra"
__email__ = "c00219805@louisiana.edu"

parser = argparse.ArgumentParser(description='Search Keys.')
parser.add_argument('--sample_name', '-sample', metavar='sample_name', \
	default='t2', help='Name of the sample on which this script should be run.')
parser.add_argument('--path', '-path', metavar='path', \
	default='//home/linc/c00219805/Research/Protien_Database//extracted_new_samples/testing/', 
	help='Directory of input sample and other files.')
parser.add_argument('--setting', '-setting', metavar='setting', \
	default='/theta29_dist35/', help='Setting with theta and distance bin values.')

def get_to_be_checked_edges(pos_arr):
	return sorted([pos_arr[0], pos_arr[1]]), sorted([pos_arr[1], pos_arr[2]]), sorted([pos_arr[0], pos_arr[2]])
def get_longest_sub_structureII(file):
	edges_vertices = {}
	checked_edges = []
	df = pd.read_csv(file, delimiter = ',', header = 0)
	df = df.sort_values(by = ['pos0','pos1','pos2'])
	print 'done reading..'
	all_pos = df[['pos0','pos1','pos2']].values.tolist()
	for first_pos in all_pos:
		to_be_checked_edges = get_to_be_checked_edges(first_pos)
		print to_be_checked_edges
		for edge in to_be_checked_edges:
			if edge not in checked_edges:
				vertices = []
				for second_pos in all_pos:
					if (edge[0] in  second_pos) and (edge[1] in second_pos):
						if (str(edge[0]) + '_' + str(edge[1]) in edges_vertices) \
							and (list(set(second_pos) - set(edge))[0] not in edges_vertices[str(edge[0]) + '_' + str(edge[1])]):
								edges_vertices[str(edge[0]) + '_' + str(edge[1])].append(list(set(second_pos) - set(edge))[0])
						else:
							edges_vertices[str(edge[0]) + '_' + str(edge[1])] = [list(set(second_pos) - set(edge))[0]]
				
				checked_edges.append(edge)
		print len(edges_vertices)


		# for second_pos in all_pos:
		# 	sorted_first_pos = sorted(first_pos)
		# 	sorted_second_pos = sorted(second_pos)
		# 	if (sorted_first_pos[0] == sorted_second_pos[0]) and (sorted_first_pos[1] == sorted_second_pos[1]):
		# 		print sorted_first_pos , sorted_second_pos
		 
		#for 
def get_path(first_pos, all_pos):
	#first_pos = [3,146,193]
		print 'first_pos: ', first_pos
		to_be_checked_edges = get_to_be_checked_edges(first_pos)
		current_triplet = first_pos
		checked_edges = []
		result = set()
		result.add((first_pos[0],first_pos[1], first_pos[2]))
		i = 0
		#got_to_outer_loop = False
		while i < len(to_be_checked_edges):	
			second_loop_count = 0		
			edge = to_be_checked_edges[i]
			if edge not in checked_edges:
				checked_edges.append(edge)
			print 'first_pos: ', first_pos,  'edge: ', edge, 'to_be_checked_edges: ', to_be_checked_edges
			for second_pos in all_pos:
				second_loop_count += 1
				got_to_outer_loop = False
				#print 'first_pos: ', first_pos,'to_be_checked_edges: ', to_be_checked_edges, 'second_pos: ', second_pos, 'current_triplet: ', current_triplet, current_triplet != second_pos, 'edge: ', edge
				if current_triplet != second_pos:
					if (edge[0] in  second_pos) and (edge[1] in second_pos):
							current_triplet = second_pos
							if (sorted([edge[1],list(set(second_pos) - set(edge))[0] ]) in checked_edges) and \
								(sorted([edge[0], list(set(second_pos) - set(edge))[0]]) in checked_edges):
									return str(first_pos[0]) + '_' + str(first_pos[1]) + '_' +str(first_pos[2]), list(result), len(list(result))
							elif sorted([edge[0], list(set(second_pos) - set(edge))[0]]) in checked_edges:
								to_be_checked_edges =  [sorted([edge[1],list(set(second_pos) - set(edge))[0] ])]
							elif sorted([edge[1],list(set(second_pos) - set(edge))[0] ]) in checked_edges:
								to_be_checked_edges =  [sorted([edge[0], list(set(second_pos) - set(edge))[0]])]	
							else:
								to_be_checked_edges =  [sorted([edge[0], list(set(second_pos) - set(edge))[0]]), sorted([edge[1],list(set(second_pos) - set(edge))[0] ])]
							result.add((edge[0], edge[1],list(set(second_pos) - set(edge))[0] ))
							
							break
			if second_loop_count == len(all_pos):
				return str(first_pos[0]) + '_' + str(first_pos[1]) + '_' +str(first_pos[2]), list(result), len(list(result))
			#i +=1

def get_longest_sub_structureIII(file):
	edges_vertices = {}
	edges_vertices_counts = {}
	
	df = pd.read_csv(file, delimiter = ',', header = 0)
	df = df.sort_values(by = ['pos0','pos1','pos2'])
	print 'done reading..'
	all_pos = df[['pos0','pos1','pos2']].values.tolist()
	#for first_pos in all_pos:
	first_pos = [2,6,13]
	triple, path, path_length = get_path(first_pos, all_pos)
	edges_vertices[triple] = path
	edges_vertices_counts[triple] = path_length
	print edges_vertices
	return pd.DataFrame(edges_vertices.items()), pd.DataFrame(edges_vertices_counts.items())

def get_longest_sub_structure(file):
	all_edges = set()
	all_vertices = set()
	df = pd.read_csv(file, delimiter = ',', header = 0)
	df = df.sort_values(by = ['pos0','pos1','pos2'])
	print 'done reading..'
	all_pos = df[['pos0','pos1','pos2']].values.tolist()
	for triple in all_pos:
		edge0, edge1, edge2 = get_to_be_checked_edges(triple)
		#print edge0, edge1, edge2
		all_edges.add(tuple(edge0))
		all_edges.add(tuple(edge1))
		all_edges.add(tuple(edge2))
		all_vertices.add(triple[0])
		all_vertices.add(triple[1])
		all_vertices.add(triple[2])
	sorted_all_vertices = sorted(list(all_vertices))
	print all_edges
	source_vertex = all_vertices[0]
	for dest_vertex in all_vertices:
		if (source_vertex,dest_vertex) in all_edges:
			print dest_vertex


if __name__ == '__main__':
	"""Executable code starts here."""
	args = parser.parse_args()

	files = glob.glob(args.path + args.sample_name + args.setting + '//*_key_search_common_keys.csv')	
	column_names = ['key', 'aa0', 'pos0', 'aa1', 'pos1', 'aa2', 'pos2', 'classT', 'theta', \
					'classL', 'distance', 'x0', 'y0', 'z0', 'x1', 'y1', 'z1', 'x2', 'y2', 'z2']

	writer = pd.ExcelWriter(args.path + args.sample_name + args.setting + \
		'common_keys/'+'//longest_substructure.xlsx', engine='xlsxwriter')
	writer_counts = pd.ExcelWriter(args.path + args.sample_name + args.setting + \
		'common_keys/'+'//longest_substructure_counts.xlsx', engine='xlsxwriter')

	for file in files:
		print ntpath.basename(file)[:4]
		get_longest_sub_structure(file)
		# path, counts = get_longest_sub_structure(file)
		# path.to_excel(writer,sheet_name=ntpath.basename(file)[:4])
		# counts.to_excel(writer,sheet_name=ntpath.basename(file)[:4])

	print('Identifying Longest Sub-Structure Completed.')

