from treetime import TreeTime
from treetime.CLI_io import create_auspice_json
from treetime.utils import parse_dates
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo

dates = parse_dates('treetime_metadata.tsv')
d = pd.read_csv('treetime_metadata.tsv', sep='\t', index_col=0)
for ri, row in d.iterrows():
    if ':' in row['treetime_date']:
        dates[ri] = [float(x) for x in row['treetime_date'].split(':')]

T = Phylo.read('tree.nwk', 'newick')
T.root_with_outgroup("BA.2")
# T.prune(T.find_any("hCoV-19/SouthAfrica/NICD-N55999/2023"))

tt = TreeTime(gtr='JC69', aln ='masked.fasta', tree=T, dates=dates, verbose=4)

long_branch = tt.tree.common_ancestor('hCoV-19/Denmark/DCGC-656489/2023', 'hCoV-19/Canada/BC-BCCDC-641728/2023')

def relative_rate(tree):
    for n in tree.find_clades():
        if n.up is None: continue
        if n==long_branch:
            n.branch_length_interpolator.gamma = 2.0
        else:
            n.branch_length_interpolator.gamma = 1.0

tt.run(infer_gtr=False, max_iter=4, branch_length_mode='joint', prune_short=False,
       resolve_polytomies=False, assign_gamma=relative_rate, clock_filter=None,
       clock_rate = 0.0005, vary_rate=0.0002, time_marginal='assign')

tt.tree.prune(tt.tree.find_any("BA.2"))
auspice = create_auspice_json(tt, timetree=True, confidence=True, seq_info=True)

with open('auspice.json', 'w') as fh:
    import json
    json.dump(auspice, fh, indent=2)

from treetime.treetime import plot_vs_years
from treetime.utils import datestring_from_numeric

def tip_labels(x):
    if x.is_terminal():
        if len(x.name)<10:
            return x.name
        else:
            return d.loc[x.name, 'country']
    else:
        return ''

def branch_labels(x):
    if x==long_branch:
        return "BA.2.86"
    else:
        return ""

fig = plt.figure(figsize=(7,5))
ax = fig.add_axes(111)
plot_vs_years(tt, confidence=[0.05, 0.95], step=1/12, label_func=tip_labels, ax=ax, branch_labels=branch_labels)

tick_labels = []
tick_positions = []
for p, x in zip(plt.xticks()[0][::2], plt.xticks()[1][::2]):
    tick_labels.append(datestring_from_numeric(float(x.get_text())))
    tick_positions.append(p)

plt.xticks(tick_positions, tick_labels, ha='right')
plt.yticks([])
plt.tick_params(rotation=30)
ax.set_axis_on()
plt.xlabel('')
plt.tight_layout()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.savefig('timetree.pdf')
