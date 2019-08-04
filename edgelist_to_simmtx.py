#!/usr/bin/python

import sys
import re

def generate_null_mtx(n):
	mtx = [] 
	for i in range(n): 
		mtx.append([0.0]*n)
	return mtx

def main():
	filename = sys.argv[1]
	n = int(sys.argv[2])
	
	sim_mtx = generate_null_mtx(n)
	
	with open(filename) as csv:
		for row in csv:
			row = re.sub(r"\n", "", row)
			triplet = row.split(',')
			i = int(triplet[0]) -1
			j = int(triplet[1]) -1
			sim_scr = float(triplet[2])
			if sim_scr > 0.0: 
				sim_mtx[i][j] = sim_scr  
				sim_mtx[j][i] = sim_scr
			else: 
				sim_mtx[i][j] = 0
				sim_mtx[j][i] = 0 
	for row in sim_mtx:
		row_str = str(row) 
		row_str = re.sub(r"[\[\],]", "", row_str) 
		print(row_str)

main()
