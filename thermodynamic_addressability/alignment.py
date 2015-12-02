import numpy as np

np.set_printoptions(linewidth=250)

def alignment(source, target):
	source = np.array(tuple(source) )
	print source
	target = np.array(tuple(target))
	previous_row = np.zeros(source.size + 1)
	print previous_row
	for s in target:
		#insertion
		current_row = previous_row + 1
		#match, mismatch
		current_row[1:] = np.minimum(current_row[1:], np.add(previous_row[:-1], source != s))
		#deletion
		current_row[1:] = np.minimum(current_row[1:], current_row[0:-1] + 1)
		previous_row = current_row
		print previous_row
	return previous_row

def local_min(values, length, tolerance):
#	local_min = np.r_[True, values[1:] < values[:-1]] & np.r_[values[:-1] <= values[1:], True] & (values <= tolerance)
	local_min =  []
	for i, v in enumerate(values):
		if i == 0:
			if v < values[i+1] and v <= tolerance:
				local_min.append(i - length)
		elif i == (len(values) - 1):
			if v <= values[i-1] and v <= tolerance:
				local_min.append(i- length)
		elif v <= values[i-1] and v < values[i+1] and v <= tolerance:
			local_min.append(i- length)
	return local_min

source = 'AAAAAAACCCCCGATTTTACACCCCCATTACATTTTTTTTTAAAAAAAA'
target = 'GATTTTACA'

row = alignment(source, target)

print 'row'

print local_min(row, len(target), 3)
