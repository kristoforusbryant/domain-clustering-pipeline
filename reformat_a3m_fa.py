#!/usr/bin/python 

import sys 
import re
import fileinput

def main():
	lines = fileinput.input()
	
	# replacing preceding deletions 
	prec_replaced = map(replacePreceding, lines) 	

	# replacing trailing deletions
	result = map(replaceTrailing, prec_replaced)
	
	# printing result 
	for t in result: 
		print(t) 

def replacePreceding(s): 
	hit = re.search("^-*", str(s)).group(0) 
	hit_len = len(hit) 
	res = '~'*hit_len + s[hit_len:]
	return res

def replaceTrailing(s):
	hit = re.search("-*$", str(s)).group(0) 
	hit_len = len(hit) 
	res = s[:len(s)-hit_len-1] + '~'*hit_len 
	return res


main()	
	
