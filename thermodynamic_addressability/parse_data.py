import json
import math
import bigfloat
import numpy
import pandas as pd
import seaborn as sns
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt



def read_data(filename, short=False):
    json_data = open(filename, 'r').read()
    raw_data = json.loads(json_data)
    seq = []
    seq_id = []
    energy = []
    for k1, v1 in raw_data.iteritems():
        if k1 == "name":
            continue
        stp_i = int(k1)
        staple = v1
        #print str(stp_i) + ' ' + v1['staple_sequence']
        for k2, v2 in staple.iteritems():
            if not k2.isdigit():
                continue
            arm_i = k2
            arm = v2
            if short and len(arm['sequence']) > 8:
				continue
            dG = arm['dG']
            min_dG = float(arm['min_dG'])
            seq.append(arm['sequence'])
            seq_id.append('stp_' + str(stp_i) + '_' + str(arm_i))
            local_min = []
            for i in range(len(dG)-1):
                if dG[i] < dG[i-1] and dG[i] <= dG[i+1]:
                    local_min.append(float(dG[i]))
            sorted_by_energy = sorted(local_min)
            energy.append(numpy.array(sorted_by_energy))
    return seq_id, seq, energy


def get_boltzmann_distribution(energy_by_arm):
    R = 8.3144621  # gas constant
    T = 293.15  # room temperature
    factor = 4184.0  # joules_per_kcal
    boltzmann_distribution = []
    for dG in energy_by_arm:
        ps = []
        total = bigfloat.BigFloat(0)
        for energy in dG:
            p = bigfloat.exp((-energy*factor)/(R*T), bigfloat.precision(1000))
            ps.append(p)
            total = bigfloat.add(total, p)
        normal_ps = []
        for p in ps:
            normal_ps.append(float(bigfloat.div(p,total)))
        boltzmann_distribution.append(numpy.array(normal_ps))
    return boltzmann_distribution
    

def plot(data):
	sns.set(style="white", palette="muted")
	f, axes = plt.subplots(2, 2, figsize=(10, 10), sharex=True)
	
	sns.despine(left=True)

	b, g, r = sns.color_palette("muted", 3)

	for i, dist in enumerate(data):
		specific = [p[0] for p in dist]
		unspecific = [sum(p[1:]) for p in dist]

		#sns.distplot(specific, kde=False, color=b, ax=axes[0])
		#sns.distplot(specific, hist=False, rug=True, color=r, ax=axes[0, 1])
		#sns.distplot(specific, hist=False, color=g, kde_kws={"shade": True}, ax=axes[1, 0])
		sns.distplot(specific, color=g, ax=axes[i,1])
		sns.distplot(unspecific, color=r, ax=axes[i,0])

	plt.setp(axes, yticks=[])
	plt.tight_layout()
	plt.savefig("/home/j3ny/experiments/energy_addressability/DeBruijn_specific_hist.png",format='png',dpi=600)
	

filename_pUC19 = 'pUC19_alpha.json'
filename_DB = 'DeBruijn_alpha.json'

_, _, energies_pUC19 = read_data(filename_pUC19)
ids, sequences, energies_pUC19_short = read_data(filename_pUC19, short=True)
_, _, energies_DB = read_data(filename_DB)
_, _, energies_DB_short = read_data(filename_DB, short=True)

pUC19_dist = get_boltzmann_distribution(energies_pUC19)
#pUC19_dist = get_boltzmann_distribution(energies_pUC19_short)
DB_dist = get_boltzmann_distribution(energies_DB)
#DB_dist = get_boltzmann_distribution(energies_DB_short)

dist = [d[0] for d in pUC19_dist]

import matplotlib
from numpy.random import randn
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(int(y))
    return s

    # The percent symbol needs escaping in latex
    #if matplotlib.rcParams['text.usetex'] == True:
    #    return s + r'$\%$'
    #else:
    #    return s + '%'

x = dist

# Make a normed histogram. It'll be multiplied by 100 later.
plt.hist(x, bins=20, normed=False)
plt.xlim(0,1)

# Create the formatter using the function to_percent. This multiplies all the
# default labels by 100, making them all percentages
formatter = FuncFormatter(to_percent)

# Set the formatter
plt.gca().yaxis.set_major_formatter(formatter)
plt.ylim(0,100)

plt.show()
