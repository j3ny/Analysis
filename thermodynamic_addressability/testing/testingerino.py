from sh import RNAfold, RNAcofold
import numpy as np

PARAMFILE = '/usr/share/ViennaRNA/dna_mathews2004.par'



def fold(seq1, seq2='', prog=RNAcofold, RNA=True):
	seq = str(seq1) if not seq2 else '%s&%s' % (seq1,seq2)
	out = prog(_in=seq, noPS=True) if RNA else prog(_in=seq, noPS=True, noconv=True, P=PARAMFILE)
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
	scaffold = ''
	with open(filename, 'r') as f:
		scaffold = f.readline().strip()
	return scaffold

def read_staples(filename):
	staples = []
	with open(filename, 'r') as f:
		for line in f:
			arms = line.strip().split('|')
			staples.append(arms) 
	return staples

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

def sliding_window_scan(arm_sequence):
	data = []
	arm_len = len(arm_sequence)
	global_align = alignment(scaffold, revcom(arm_sequence))
	for i in range(100):#len(scaffold) - arm_len):
		target = scaffold[i:i+arm_len]
		levenstein = alignment(revcom(arm_sequence), target)[-1]
		g_align = global_align[i+arm_len]
		structure, energy = fold(arm_sequence, target, RNA=False)
		#print i, target, levenstein, g_align, structure, energy
		data.append((target, levenstein, g_align, structure, energy))
	return data


def local_min(values, length, tolerance):
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


path = ''
scaffold_filename = path + 'scaffold.txt'
staples_filename = path + 'DXstaples.txt'
scaffold = read_scaffold(scaffold_filename)

def main():
	staples = read_staples(staples_filename)

	staples = staples[:1]

	staple_data = {}
	for i, staple in enumerate(staples):
		for arm in staple:
			print str(i+1) + ' out of ' + str(len(staples))
			staple_data[i] = sliding_window_scan(arm)

	print 'done analysing'

	raw_filename = 'data_raw.csv'
	with open(raw_filename, 'w') as log:
		for i, staple in enumerate(staples):
			arm = staple[0]
			log.write(arm + '\n')
			for data in staple_data[i]:
				#print data
				log.write(str(data) + '\n')

	filename = 'data_local_min.csv'
	with open(filename, 'w') as log:
		for i, staple in enumerate(staples):
			arm = staple[0]
			log.write(arm + '\n')
			log.write('target,levenstein,global_align,structure,energy,energy_minimum_right,energy_minimum_left,alignment_minimum\n')
			align = alignment(scaffold, revcom(arm))
			loc_min = local_min(align, len(arm), tolerance=8)
			log.write(str(align[:100]) + '\n')
			log.write(str(loc_min[:100]) + '\n')
			for j in range(len(staple_data[i])):
				t, l, glob, s, energy = staple_data[i][j]
				e_prev = staple_data[i][j-1][4] if j > 0 else 0.0
				e_next = staple_data[i][j+1][4] if j < len(staple_data[i]) - 1 else 0.0
				e_minimum_r = '1' if (energy < e_prev and energy <= e_next) else '0'
				e_minimum_l = '1' if (energy <= e_prev and energy < e_next) else '0'
				g_minimum = '1' if (j in loc_min) else '0'
				log.write('{},{},{},{},{},{},{},{}\n'.format(t, int(l), int(glob), s, energy, e_minimum_r, e_minimum_l, g_minimum))

if __name__ == "__main__":
    main()







