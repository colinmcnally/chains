# print filename base, for moving files to TMP_DIR
import chainmodels
import sys


runindex = int(sys.argv[1])-1 #job array indexes start at 1

keys = list(chainmodels.runs_ext2.keys())
keys = sorted(keys)
hashkey = keys[runindex]
print(chainmodels.filenames_ext2[hashkey][:-3])
