import sys
import numpy as np


in_filename = 'S.txt.filtered'

argc = len(sys.argv)
argcc = 1
while argcc < argc:
    if sys.argv[argcc] == '-i':
        argcc += 1
        in_filename = sys.argv[argcc]
    argcc += 1

info = []
with open(in_filename) as f:
    done = False
    while not done:
        line = f.readline()
        if line == '':
            done = True
        else:
            if line[0] == '@':
                tokens = line.rstrip().split(' ')
                qscore = float(tokens[3])
                pos = int(tokens[5])
                read = f.readline()
                match = f.readline().rstrip()
                perfect = match.find(' ')
                if perfect == -1:
                    perfect = len(match)
                info.append([qscore, pos, perfect])
info = np.array(info)

# analysis

for i in range(5,14):
    print('len %d perfect: %d' % (i, np.count_nonzero(info[:,2] == i)))

