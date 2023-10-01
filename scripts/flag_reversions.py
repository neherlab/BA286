"""
Script to flag reversions to reference in the output Auspice JSON
Takes:
- Reference fasta
- Auspice JSON
Outputs:
- TSV with node/strain and reversion
"""
#%%
from Bio import SeqIO
from Bio.Seq import Seq
import json

#%%

# Load reference as string
ref = str(SeqIO.read("config/reference.fasta", "fasta").seq)
# print(ref)
# %%
# Get .tree from auspice json
tree=json.load(open("auspice/BA.2.86.json"))["tree"]
# print(tree)
# %%
# Tree is recursive with children under .children
# Mutations are a list under .branch_attrs.mutations.nuc
# Reversion is defined as a mutation to the reference state
# 

def get_reversions(tree, ref):
    reversions = []
    for node in tree:
        if "children" in node:
            reversions += get_reversions(node["children"], ref)
        node_attrs = node.get("node_attrs", {})
        num_date = node_attrs.get("num_date", {}).get("value", 0)
        branch_attrs = node.get("branch_attrs", {})
        mutations = branch_attrs.get("mutations", {})
        nucs = mutations.get("nuc", [])
        for nuc in nucs:
            pre, pos, post = nuc[0], int(nuc[1:-1]), nuc[-1]
            if post == ref[pos-1]:
                # Recent samples with node_attrs.num_date.value >=2023
                if num_date >= 2023:
                    reversions.append((node["name"], nuc))
    return reversions

reversions = get_reversions(tree["children"], ref)

# %%

for node, reversion in reversions:
    print(f"{node}\t{reversion[1:-1]}")
# %%
