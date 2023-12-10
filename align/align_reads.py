import sys
import numpy as np
import math
import editdistance
import matplotlib.pyplot as plt

verbose = 0

# simple function to load a fasta file
def load_fasta(filename):
    genome = ''
    description = ''
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.rstrip()
            if line[0] == '>':
                description = line[1:]
            else:
                genome += line

    return genome, description

def load_fastq(filename):
    reads = []
    with open(filename) as f:
        done = False
        while not done:
            line = f.readline()
            if line == '':
                done = True
            elif line[0] == '@':
                reads.append(f.readline().rstrip())
    return reads

def reverse_comp(read):
    rcomp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    rread = ''
    for c in read[::-1]:
        rread += rcomp[c]
    return rread

def bases2vals(bases):
    valdict = {'A':'0', 'C':'1', 'G':'2', 'T':'3'}
    vals = ''
    for base in bases:
        vals += valdict[base]
    return vals

# generates an array of size 4^hash_len and each entry contains a list of genome positions matching the hash
def generate_hashlist(genome, hash_len):
    hash_list = []
    genome_len = len(genome)
    num_hashs = int(math.pow(4,hash_len))

    for i in range(num_hashs):
        hash_list.append([])

    genome_position = 0
    while genome_position < (genome_len-hash_len):
        if pos_list and genome_position not in pos_list:
            genome_position += 1
            continue
        index = int(bases2vals(genome[genome_position:genome_position+hash_len]), 4)
        hash_list[index].append(genome_position)
        genome_position += 1

    return hash_list

def map_read(genome, read, hash_list, hash_len):
    # for each genome position in the hash list, compute edit distance for read, return the position & edit distance of the best match
    bestdist = -1
    bestpos = -1
    readlen = len(read)
    hash_index = int(bases2vals(read[:hash_len]), 4) # cool python way to convert a base-4 string to a base-10 integer
    for pos in hash_list[hash_index]:
        dist = editdistance.eval(genome[pos:pos+readlen], read)
        if bestdist == -1 or dist < bestdist:
            bestdist = dist
            bestpos = pos

    return bestpos,bestdist

def qualscore(readlen, num_errors):
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

    return (bestq,best_q_len,perfect_len)

def top_n_by_len(N, L):
    forced_info = []
    for i in range(len(info)):
        # start calculating the Q-Score at length L (our minimum) and continue to the end of the read to get the best Q-Score
        forced_info.append(scoremin(info[i][0], info[i][1], L))
    forced_info = np.array(forced_info)

    top_by_q_i = np.flip(np.argsort(forced_info[:,0]))
    top_by_q = forced_info[top_by_q_i]
    print('%d %dQ%.2f' % (N, L, np.mean(top_by_q[:N,0])))
    return top_by_q_i

# set some defaults
ref_name = 'phix174.fasta'
in_filename = None
out_filename = None
hash_len = 4
score_cutoff = 0.0
score_cutofflen = 0
direction_filter = 0
pos_list = None
want_plots = False

# parse cmd-line args
argcc = 1
argc = len(sys.argv)
while argcc < argc:
    if sys.argv[argcc] == '--ref':
        argcc += 1
        ref_name = sys.argv[argcc]
    if sys.argv[argcc] == '--in':
        argcc += 1
        in_filename = sys.argv[argcc]
    if sys.argv[argcc] == '--out':
        argcc += 1
        out_filename = sys.argv[argcc]
    if sys.argv[argcc] == '--hl':
        argcc += 1
        hash_len = int(sys.argv[argcc])
    if sys.argv[argcc] == '--filter':
        argcc += 1
        score_cutoff = float(sys.argv[argcc])
        argcc += 1
        score_cutofflen = int(sys.argv[argcc])
    if sys.argv[argcc] == '--direction':
        argcc += 1
        direction_filter = int(sys.argv[argcc])
    if sys.argv[argcc] == '--poslist':
        argcc += 1
        pos_list = [int(i) for i in sys.argv[argcc].split(',')]
    if sys.argv[argcc] == '--plots':
        want_plots = True
    if sys.argv[argcc] == '-v':
        verbose += 1
    argcc += 1

# load up our reads
if in_filename is None:
    print('missing --in readfile?')
    exit(0)
reads = load_fastq(in_filename)
print('loaded %d reads' % len(reads))

# load up the fasta file
ref,description = load_fasta(ref_name)
if verbose > 0:
    print('loaded ref: %s\nlength: %d\n' % (description, len(ref)))

if verbose > 1:
    print('first few reads:')
    for r in range(20):
        print('%s' % reads[r])

# generate hash list from genome and map reads
print('generating hash list...')
hash_list = generate_hashlist(ref, hash_len)
ref_r = reverse_comp(ref)
hash_list_r = generate_hashlist(ref_r, hash_len)

if verbose > 2:
    print('reference:\n%s' % ref)
    print('rcomp ref:\n%s' % ref_r)
    print('reference hash list:\n%s' % str(hash_list))
    print('rcomp ref hash list:\n%s' % str(hash_list_r))

outfile = None
if out_filename:
    outfile = open(out_filename, 'w') 

# stores reads & ref so we can track some mapping stats
info = []

print('mapping reads...')
for i, read in enumerate(reads):
    rcomp = False
    # we map to both the forward and the reverse complement of the read, to see which is the best match
    pos,dist = map_read(ref, read, hash_list, hash_len)
    pos_r,dist_r = map_read(ref_r, read, hash_list_r, hash_len)

    if pos > -1 or pos_r > -1: # if it mapped
        if pos > -1 and pos_r > -1: # if it mapped both directions, pick the best one
            if dist_r < dist:
                pos = pos_r
                rcomp = True
        elif pos_r > -1: # else if it just mapped to the rcomp, use that one
            pos = pos_r
            rcomp = True
                         # else just stick with the forward read pos & dist

        ref_read = ref_r[pos:pos+len(read)] if rcomp else ref[pos:pos+len(read)]
        if len(ref_read) == len(read):
            outtxt = 'read: %s  ref: %s  rcomp: %s  pos: %d  dist: %d' % (read, ref_read, rcomp, pos, dist)
            if outfile:
                outfile.write('%s\n' % outtxt)
            else:
                print(outtxt)
            read_name = 'read_' + str(i)

            # store reads & ref so we can calc some stats later
            info.append((read, ref_read, rcomp, pos, i))
        else:
            if verbose > 0:
                print('failed to align read: %s to valid reference position' % read)

if outfile:
    outfile.close()

print('mapped %d out of %d reads' % (len(info), len(reads)))

#####################
# evaluations
#####################

# calculate and plot coverage
# also use this to filter out junk
info_filtered = []
cov = np.zeros(len(ref))
cov_filtered = np.zeros(len(ref))
cov_starts = np.zeros(len(ref))
num_forward = 0
num_rcomp = 0
num_filtered = 0
for read in info:
    score = scoremin(read[0], read[1], score_cutofflen)
    if score[0] >= score_cutoff:
        start = read[3]
        if read[2] is True: # rcomp:
            start -= len(read[0])
            num_rcomp += 1
            direction = -1
        else:
            num_forward += 1
            direction = 1

        if direction_filter == 0 or direction_filter == direction:
            info_filtered.append(read)
            cov_filtered[start:start+len(read[0])] += 1
    cov[start:start+len(read[0])] += 1
    cov_starts[start] += 1

num_filtered = len(info_filtered)
print('num filtered: %d  forward: %d  rcomp: %d' % (num_filtered, num_forward, num_rcomp))

ymax = np.max(cov)
ymax = round((ymax+99)/100) * 100
if ymax < 100:
    ymax = 100

plt.figure('coverage')
plt.plot(cov)
ax = plt.gca()
ax.set_ylim([0, ymax])
plt.savefig('coverage.png')

plt.figure('starts')
plt.plot(cov_starts)
plt.savefig('starts.png')

plt.figure('filtered coverage')
plt.plot(cov_filtered)
ax = plt.gca()
ax.set_ylim([0, ymax])
plt.savefig('filtered_coverage.png')


print('filtered')
info = info_filtered
top_n_by_len(50 if num_filtered > 50 else num_filtered, 6)
top_n_by_len(100 if num_filtered > 100 else num_filtered, 6)
top_by_q_i = top_n_by_len(50 if num_filtered > 50 else num_filtered, 10)
top_n_by_len(100 if num_filtered > 100 else num_filtered, 10)
top_n_by_len(100 if num_filtered > 100 else num_filtered, 12)

cov_top_n = np.zeros(len(ref))
for i in range(50):
    read = info[top_by_q_i[i]]
    start = read[3]
    if read[2] is True: # rcomp:
        start -= len(read[0])
        num_rcomp += 1
        direction = -1
    else:
        num_forward += 1
        direction = 1

    if direction_filter == 0 or direction_filter == direction:
        cov_top_n[start:start+len(read[0])] += 1

plt.figure('coverage of top 50Q10')
plt.plot(cov_top_n)
plt.savefig('coverage_top_50Q10.png')


counts= {}
counts['A'] = 0
counts['C'] = 0
counts['G'] = 0
counts['T'] = 0

avg_scores = 0.0
avg_scores5 = 0.0
avg_scores10 = 0.0
filtered_filename = out_filename + '.filtered'
with open(filtered_filename, 'w') as f:
    for read in info_filtered:
        bars = ''
        for i in range(len(read[0])):
            if read[0][i] == read[1][i]:
                bars += '|'
            else:
                bars += ' '
            counts[read[0][i]] += 1
        f.write('@read: %d q-score: %.2f pos: %d rcomp: %s\n%s\n%s\n%s\n' % (read[4], scoremin(read[0], read[1], 13)[0], read[3], read[2], read[0], bars, read[1]))
        avg_scores += scoremin(read[0], read[1], 13)[0]
        avg_scores5 += scoremin(read[0], read[1], 5)[0]
        avg_scores10 += scoremin(read[0], read[1], 10)[0]

avg_scores /= len(info_filtered)
avg_scores5 /= len(info_filtered)
avg_scores10 /= len(info_filtered)
print('%d filtered reads with avg q-score at full length: %.2f  len5: %.2f  len10: %.2f' % (len(info_filtered), avg_scores, avg_scores5, avg_scores10))

print('base counts: A: %d  C: %d  G: %d  T: %d\n' % (counts['A'], counts['C'], counts['G'], counts['T']))

if want_plots:
    plt.show()

