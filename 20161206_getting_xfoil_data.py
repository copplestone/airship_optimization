import glob
import csv
from subprocess import *
from itertools import product
from os import remove
from multiprocessing import Pool


# Not the first instance as this opens the program itself whereas we want to access an executable

cmd = '/Applications/Xfoil.app/Contents/Resources/xfoil'
xfoil = Popen(cmd, stdin=PIPE, stdout=PIPE, shell=True, bufsize=-1, universal_newlines=True)
xfoil.stdin.write("\n".join(["NACA 0012","OPER","alfa 0","CPWR","test_1240"]))
xfoil.stdin.close()
print "Data stored"

with open('test_1240', 'rb') as csvfile:
	CPx = csv.reader(csvfile, delimiter=' ', quotechar='|')
	data = []
	for row in CPx:
		new_info = row[5],row[-1]
		data.append(new_info)
		print new_info

