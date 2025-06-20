import sys

with open(sys.argv[1],"r") as input:
	for line in input:
		line = line.split("\t")[0] + "\t" + line.split("\t")[1].split('.')[0] 
		print(line)
