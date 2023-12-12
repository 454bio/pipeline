import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import math

def qualscore(readlen, num_errors):
    if num_errors == 0:
        return 30.0
    qual_score = -10.0 * math.log10(num_errors/readlen)
    return qual_score

def scoremin(seq, ref, minlen):
    bestq = 0
    best_q_len = 0
    perfect_len = 0

    errors = 0
    ls = len(seq)
    lr = len(ref)
    if lr < ls:
        ls = lr

    q = 0

    for i in range(minlen):
        if seq[i] != ref[i]:
            errors += 1
        if perfect_len == 0 and errors == 1:
            perfect_len = i
    if errors > 0:
        q = qualscore(minlen, errors)
        bestq = q
        best_q_len = minlen
    else:
        perfect_len = minlen
        bestq = 30.0

    for i in range(minlen,ls):
        if seq[i] != ref[i]:
            errors += 1
        if errors > 0:
            q = qualscore(i+1, errors)
        else:
            q = 0
        if perfect_len == 0 and errors == 1:
            perfect_len = i
        if q > bestq:
            bestq = q
            best_q_len = i+1

    #return (bestq,best_q_len,perfect_len)
    return bestq

def lenq(seq, ref, minq):
    ls = len(seq)
    lr = len(ref)
    if lr < ls:
        ls = lr
    bestLen = 0
    errors = 0
    for i in range(ls):
        if seq[i] != ref[i]:
            errors += 1
        qscore = qualscore(i+1, errors)
        if qscore >= 10.0:
            bestLen = i+1
    return bestLen


in_filename = 'S.txt.filtered'
show_plots = False

argc = len(sys.argv)
argcc = 1
while argcc < argc:
    if sys.argv[argcc] == '-i':
        argcc += 1
        in_filename = sys.argv[argcc]
    if sys.argv[argcc] == '--plots':
        show_plots = True
    argcc += 1

run_name = in_filename.split('.')[0]

# load up the "filtered" read info

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
                spot_id = int(tokens[1])
                qscore = float(tokens[3])
                pos = int(tokens[5])
                read = f.readline()
                readlen = len(read)
                match = f.readline().rstrip()
                perfect = match.find(' ')
                if perfect == -1:
                    perfect = len(match)
                ref = f.readline().rstrip()
                q10 = scoremin(read, ref, 10)
                lenq10 = lenq(read, ref, 10)
                info.append([spot_id, qscore, q10, pos, perfect, readlen, lenq10])

# open the fastq file, and for each read we can find in the filtered info, include additional metrics
with open(run_name + ".fastq") as f:
    done = False
    info_index = 0
    while not done:
        line = f.readline()
        if line == '':
            done = True
        else:
            if line[0] == '@':
                tokens = line.rstrip().split(' ')
                if len(tokens) > 1:
                    spot_id = int(tokens[0].split('_')[1])
                    if spot_id == info[info_index][0]:
                        j = json.loads(tokens[1])
                        info[info_index].append(j['err'])
                        _ = f.readline()
                        _ = f.readline()
                        qual_text = f.readline()
                        # append the first few quality scores
                        for i in range(10):
                            qval = ord(qual_text[i]) - 33
                            info[info_index].append(qval)
                        info_index += 1
                        if info_index >= len(info):
                            done = True

# info contains spot_id, qscore, Q10 score, pos of map, perfect len, readlen, lenq10, phase fit err, 8: pos 0 err, 9: pos 1 err, ..., pos 9 err
info = np.array(info)

# analysis

readlen = int(info[0,5]) # just pull length from first read and assume they are all this way - hack for now
num = len(info)

for i in range(5,readlen+1):
    print('len %d perfect: %d' % (i, np.count_nonzero(info[:,4] == i)))

# try and pull out which template we mapped to, based on position
max_pos = int(np.max(info[:,3]))
pos_array = np.zeros(max_pos+1);
for i in range(num):
    pos_array[int(info[i,3])] += 1
for i in range(max_pos+1):
    if pos_array[i] > 0:
        print('pos %d: %d reads' % (i, int(pos_array[i])))

if info.shape[1] > 7:
    plt.figure('perfect vs qual')
    plt.scatter(info[:,4], info[:,7])

'''
plt.figure('quals per position')
for i in range(10):
    plt.scatter([i]*num, info[:,i+7])
'''

if info.shape[1] > 7:
    fix,axes = plt.subplots(nrows=1, ncols=10)
    for i in range(10):
        axes[i].boxplot(info[:,i+8], vert=True, patch_artist=True, labels=['pos ' + str(i+1)])
        axes[i].set_ylim(0, 30)

plt.figure('10Q10 all')
plt.hist(info[:,2], bins=21)
print('all: %d 10Q%.2f' % (num, np.mean(info[:,2])))

plt.figure('len >= Q10')
plt.hist(info[:,6], bins=range(16))

filtered = []
for i in range(num):
    if info[i,2] >= 10:
        filtered.append(info[i]) # wow, such a hack - I need a numpy way to select, which I know exists
filtered = np.array(filtered)
plt.figure('cheat 10Q10 filtered')
plt.hist(filtered[:,2])
print('cheat filtered: %d 10Q%.2f' % (len(filtered), np.mean(filtered[:,2])))

whitelist_name = run_name + '.whitelist.txt'
with open(whitelist_name, 'w') as f:
    for i in range(len(filtered)):
        val = filtered[i][0] # white-list needs to be zero-based for the C++ Basecaller
        f.write('%d\n' % val)
    

if info.shape[1] > 7:
    #filter_params=[12,12,12,1]
    #filter_params=[9,5,9,9]
    #filter_params=[9,9,9,9]
    filter_params=[10,10,10,10]
    real_filtered = []
    for i in range(num):
        if info[i,8] >= filter_params[0] and info[i,9] >= filter_params[1] and info[i,10] >= filter_params[2] and info[i,11] >= filter_params[3]:
            real_filtered.append(info[i])
    real_filtered = np.array(real_filtered)
    plt.figure('real 10Q10 filtered')
    plt.hist(real_filtered[:,2], bins=21)
    print('real filtered: %d 10Q%.2f' % (len(real_filtered), np.mean(real_filtered[:,2])))

if show_plots:
    plt.show()

