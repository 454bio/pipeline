import sys
import numpy as np
import json
import matplotlib.pyplot as plt

filename = 'reads.fastq'
err_cutoff = 0

argcc = 1
argc = len(sys.argv)
while argcc < argc:
    if sys.argv[argcc] == '-i':
        argcc += 1
        filename = sys.argv[argcc]
    if sys.argv[argcc] == '--cutoff':
        argcc += 1
        err_cutoff = float(sys.argv[argcc])
    argcc += 1

# read in the optional json params from the info line of each read, accumulate stats
vals = []
with open(filename) as f:
    lines = f.readlines()
    for line in lines:
        if line[0] == '@':
            tokens = line.strip().split(' ')
            if len(tokens) >= 2:
                stats = json.loads(tokens[1])
                if 'ie' in stats:
                    vals.append([stats['ie'], stats['cf'], stats['dr'], stats['err']])
                else:
                    vals.append(stats['err'])

vals = np.array(vals)
print('shape: %s' % str(vals.shape))

if vals.ndim == 1: # we just have err, show avg
    errvals = vals
else:
    errvals = vals[:,3]

print('error mean: %f median: %f' % (np.mean(errvals), np.median(errvals)))

if vals.ndim  > 1:
    if err_cutoff == 0:
        err_cutoff = np.median(errvals) * 2
    print('using error cutoff value: %f' % err_cutoff)
    ievals = vals[vals[:,3] < err_cutoff,0]
    cfvals = vals[vals[:,3] < err_cutoff,1]
    drvals = vals[vals[:,3] < err_cutoff,2]
    print('ie mean: %f median: %f' % (np.mean(ievals), np.median(ievals)))
    print('cf mean: %f median: %f' % (np.mean(cfvals), np.median(cfvals)))
    print('dr mean: %f median: %f' % (np.mean(drvals), np.median(drvals)))


plt.figure('err')
plt.hist(errvals.clip(min=0.0, max=1.0), bins=51, range=[0.0, 1.0])

plt.figure('err_all')
plt.hist(errvals, bins=101, range=[0.0, 10.0])

plt.figure('ie')
plt.hist(ievals, bins=21)

plt.figure('cf')
plt.hist(cfvals, bins=21)

plt.figure('dr')
plt.hist(drvals, bins=21)

plt.show()

