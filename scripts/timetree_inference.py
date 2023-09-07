import argparse
import json
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo
from treetime import TreeTime
from treetime.CLI_io import create_auspice_json
from treetime.utils import parse_dates



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--metadata', default='builds/metadata.tsv')
    parser.add_argument('--tree', default='builds/tree.nwk')
    parser.add_argument('--output-tree')
    parser.add_argument('--output-node-data')
    parser.add_argument('--alignment', default='builds/aligned.fasta')
    parser.add_argument('--figure', default='figures/timetree.pdf')
    args = parser.parse_args()

    dates = parse_dates(args.metadata)
    d = pd.read_csv(args.metadata, sep='\t', index_col=0)
    d.region[d.region=='Africa'] = 'Other Africa'
    d.region[(d.country=='South Africa') | (d.country=='Eswatini') | (d.country=='Namibia')] = 'Southern Africa'

    # for ri, row in d.iterrows():
    #     if ':' in row['treetime_date']:
    #         dates[ri] = [float(x) for x in row['treetime_date'].split(':')]

    T = Phylo.read(args.tree, 'newick')
    T.root_with_outgroup("BA.2")
    # T.prune(T.find_any("hCoV-19/SouthAfrica/NICD-N55999/2023"))

    tt = TreeTime(gtr='JC69', aln =args.alignment, tree=T, dates=dates, verbose=4)

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

    if args.output_tree:
        Phylo.write(tt.tree, args.output_tree, 'newick')

    node_data = {"nodes":{}}
    for n in tt.tree.find_clades():
        conf = tt.get_max_posterior_region(n, fraction=0.9)
        node_data['nodes'][n.name] = {'branch_length': len([p for a,p,d in n.mutations if a in 'ACGT' and d in 'ACGT']),
                                      'clock_length': n.clock_length,
                                      'numdate':n.numdate,
                                      'date': n.date,
                                      'num_date_confidence':[float(conf[0]), float(conf[1])]}

    if args.output_node_data:
        import json
        with open(args.output_node_data, 'w') as fh:
            json.dump(node_data, fh, indent=1)

    from treetime.treetime import plot_vs_years
    from treetime.utils import datestring_from_numeric

    def tip_labels(x):
        return ''

    def branch_labels(x):
        if x==long_branch:
            return "BA.2.86"
        else:
            return ""

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes(111)
    fs=12
    from matplotlib.colors import to_hex
    region_colors = {r:to_hex(f'C{i}') for i,r in enumerate(d.region.unique())}

    for y,n in enumerate(tt.tree.get_terminals()):
        n.y = y
        n.color='#888888'
        #n.color = region_colors[d.loc[n.name, 'region']]

    for n in tt.tree.get_nonterminals():
        n.color='#888888'

    plot_vs_years(tt, confidence=[0.05, 0.95], step=1/12, label_func=tip_labels, ax=ax,
                  branch_labels=branch_labels, selective_confidence= lambda x: x==long_branch)

    tt.tree.root.x = tt.tree.root.branch_length
    for n in tt.tree.get_nonterminals(order='preorder'):
        for c in n:
            c.x = n.x + c.branch_length


    for region in sorted(region_colors.keys()):
        if region=='?': continue
        tips = [n for n in tt.tree.get_terminals() if d.loc[n.name, 'region']==region]
        plt.scatter([n.x for n in tips],  [n.y for n in tips],
                    c=region_colors[region], label=region, zorder=2, s=20)

    plt.legend(fontsize=fs)
    tick_labels = []
    tick_positions = []
    for p, x in zip(plt.xticks()[0][::2], plt.xticks()[1][::2]):
        tick_labels.append(datestring_from_numeric(float(x.get_text())))
        tick_positions.append(p)

    plt.xticks(tick_positions, tick_labels, ha='right')
    plt.yticks([])
    plt.tick_params(rotation=30, labelsize=fs)
    ax.set_axis_on()
    plt.xlabel('')
    plt.tight_layout()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.savefig(args.figure)
