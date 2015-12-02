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
    
def alignment(source, target):
	source = np.array(tuple(source) )
	target = np.array(tuple(target))
	previous_row = np.zeros(source.size + 1)
	for s in target:
		#insertion
		current_row = previous_row + 1
		#match, mismatch
		current_row[1:] = np.minimum(current_row[1:], np.add(previous_row[:-1], source != s))
		#deletion
		current_row[1:] = np.minimum(current_row[1:], current_row[0:-1] + 1)
		previous_row = current_row
	return previous_row

def local_min(values, length, tolerance):
#	local_min = np.r_[True, values[1:] < values[:-1]] & np.r_[values[:-1] <= values[1:], True] & (values <= tolerance)
	local_min =  []
	for i, v in enumerate(values):
		if i == 0:
			if v < values[i+1] and v <= tolerance:
				local_min.append(0)
		elif i == (len(values) - 1):
			if v <= values[i-1] and v <= tolerance:
				local_min.append(i - length)
		elif v <= values[i-1] and v < values[i+1] and v <= tolerance:
			local_min.append(i - length) if i - length >= 0 else local_min.append(0)
	return local_min

def scan_arm(arm_sequence):
	target = revcom(arm_sequence)
	align = alignment(scaffold, target)
	loc_min = local_min(align, len(arm_sequence), tolerance=8)
	arm_data = []
	for i in loc_min:
		region = scaffold[i:i+len(arm_sequence)]
		if region == '':
			print arm_sequence, target, i
		structure, energy = fold(arm_sequence, region, RNA=False)#[1]
		#lev = levenshtein_dist(revcom(arm_sequence), target)
		arm_data.append((energy, region, structure))
	sorted_by_energy = sorted(arm_data, key=lambda tup: tup[0])
	return sorted_by_energy
    
def analyze(scaffold, staples, filename):
    with open(filename, 'w') as output:  
        for i, staple in enumerate(staples):
            print i
            data = scan_staple(staple)
            for arm, arm_data in data:
                output.write(arm + '\n')
                for energy, target, structure in arm_data:
                    output.write("{},{},{}\n".format(energy, target, structure))

path = 'data/raw/sequences/M13_medium/'
scaffold_filename = path + 'scaffold.txt'
staples_filename = path + 'DXstaples.txt'

scaffold = read_scaffold(scaffold_filename)
staples = read_staples(staples_filename)

analyze(scaffold, staples, 'data/verification/M13_medium_tol8.csv')

path = 'data/raw/sequences/DBS_medium/'
scaffold_filename = path + 'scaffold.txt'
staples_filename = path + 'DXstaples.txt'

scaffold = read_scaffold(scaffold_filename)
staples = read_staples(staples_filename)

analyze(scaffold, staples, 'data/verification/DBS_medium_tol8.csv')

