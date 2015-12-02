import numpy as np

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
	
print alignment('TTCTGTGAACTG', 'TTCTGTGA')
