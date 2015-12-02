import numpy as np
from collections import OrderedDict
from sh import RNAfold, RNAcofold

PARAMFILE = '/usr/share/ViennaRNA/dna_mathews2004.par'

def fold(seq1, seq2='', prog=RNAcofold, RNA=True):
    seq = str(seq1) if not seq2 else '%s&%s' % (seq1,seq2)
    out = prog(_in=seq, noPS=True) if RNA else prog(_in=seq, noPS=True, noconv=True, dangles=0, P=PARAMFILE)
    seq, fold, _ = out.split('\n')
    structure = str(fold.split()[0])
    energy = float(fold.partition(' ')[2][1:-1].strip())
    return structure, energy

def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)
def revcom(s):
    return complement(s[::-1])

def read_scaffold(filename):
    f = open(filename, 'r')
    return f.readline().strip()

def read_staples(filename):
    f = open(filename, 'r')
    staples = []
    for line in f:
        arms = line.strip().split('|')
        staples.append(arms) 
    return staples

def scan_staple(staple):
    data = []
    for arm in staple:
        data.append((arm, scan_arm(arm)))
    return data



def scan_arm(arm_sequence):  
    #sliding window over scaffold
    levenshtein = []
    for i in range(len(scaffold) - len(arm_sequence) + 1):
        target = scaffold[i:i+len(arm_sequence)]
        lev = levenshtein_dist(revcom(arm_sequence), target)
        levenshtein.append(lev)
    lev_loc_min = []
    threshold = 3
    for i in range(len(levenshtein)):
        l_prev = levenshtein[i-1] if i-1>0 else len(arm_sequence)
        l = levenshtein[i]
        l_next = levenshtein[i+1] if i+1<len(levenshtein) else len(arm_sequence)
        if l < l_prev and l <= l_next and l <= threshold:
            lev_loc_min.append(i)
    arm_data = []
    for i in lev_loc_min:
        target = scaffold[i:i+len(arm_sequence)]
        structure, energy = fold(arm_sequence, target, RNA=False)#[1]
        hamming = hamming_distance(revcom(arm_sequence), target)
        lev = levenshtein_dist(revcom(arm_sequence), target)
        arm_data.append((energy, hamming, lev, target, structure))
    #consider only local minima sites
    #arm_data = []
    #local_min = []
    #for i in range(len(arm_data)):
    #    e_prev = arm_data[i-1][0] if i-1>0 else 0.0
    #    e = arm_data[i][0]
    #    e_next = arm_data[i+1][0] if i+1<len(arm_data) else 0.0
    #    if e < e_prev and e <= e_next:
    #        local_min.append(arm_data[i])
    sorted_by_energy = sorted(arm_data, key=lambda tup: tup[0])
    return sorted_by_energy

def levenshtein_dist(source, target):
    if len(source) < len(target):
        return levenshtein_dist(target, source)
    # So now we have len(source) >= len(target).
    if len(target) == 0:
        return len(source)
    # We call tuple() to force strings to be used as sequences
    # ('c', 'a', 't', 's') - numpy uses them as values by default.
    source = np.array(tuple(source))
    target = np.array(tuple(target))
    # We use a dynamic programming algorithm, but with the
    # added optimization that we only need the last two rows
    # of the matrix.
    previous_row = np.arange(target.size + 1)
    for s in source:
        # Insertion (target grows longer than source):
        current_row = previous_row + 1
        # Substitution or matching:
        # Target and source items are aligned, and either
        # are different (cost of 1), or are the same (cost of 0).
        current_row[1:] = np.minimum(
                current_row[1:],
                np.add(previous_row[:-1], target != s))
        # Deletion (target grows shorter than source):
        current_row[1:] = np.minimum(
                current_row[1:],
                current_row[0:-1] + 1)
        previous_row = current_row
    return previous_row#[-1]

def hamming_distance(source, target):
    assert len(source) == len(target)
    count = 0
    for p, t in zip(source, target):
        count += 1 if p != t else 0
    return count
    
def analyze(scaffold, staples, filename):
    with open(filename, 'w') as output:  
        for i, staple in enumerate(staples):
            print i
            data = scan_staple(staple)
            for arm, arm_data in data:
                output.write(arm + '\n')
                for energy, hamming, levenshtein, target, structure in arm_data:
                    output.write("{},{},{},{},{}\n".format(energy, hamming, levenshtein, target, structure))

#path = 'data/raw/sequences/M13_medium/'
#scaffold_filename = path + 'scaffold.txt'
#staples_filename = path + 'DXstaples.txt'

#scaffold = read_scaffold(scaffold_filename)
#staples = read_staples(staples_filename)

#analyze(scaffold, staples, 'data/verification/M13_medium_raw.csv')

#path = 'data/raw/sequences/DBS_medium/'
#scaffold_filename = path + 'scaffold.txt'
#staples_filename = path + 'DXstaples.txt'

#scaffold = read_scaffold(scaffold_filename)
#staples = read_staples(staples_filename)

#analyze(scaffold, staples, 'data/verification/DBS_medium_raw.csv')

print levenshtein_dist('12345678901234567890', '1234567890')
