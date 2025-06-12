"""
Post-process auspice.json for publication svgs
## All
- Rotate branch with C9866T to the top

## Basal inset
- Prune all children from branch with C897A
- Set divergence of this branch to 11

"""
#%%
import json

#%%
auspice=json.load(open('builds/BA.3.2.json'))
modified = json.load(open('builds/BA.3.2.json'))

#%%
# Locate node that has C9866T in nuc
def get_node_with_mutation(node, mutation, max_div=2):
    if 'branch_attrs' in node:
        branch_attrs = node.get('branch_attrs', {})
        mutations = branch_attrs.get('mutations', {})
        nuc = mutations.get('nuc', [])
        if any(mutation == m for m in nuc) and node["node_attrs"]["div"] < max_div:
            return {"child": node}
    if 'children' in node:
        for child in node['children']:
            node_with_mutation = get_node_with_mutation(child, mutation, max_div=max_div)
            if node_with_mutation is None:
                continue
            if "parent" in node_with_mutation:
                return node_with_mutation
            elif "child" in node_with_mutation:
                return {"parent": node, "child": node_with_mutation["child"]}
    return None

# Recursively prune child
def reorder_child(parent, child_name):
    if 'children' in parent:
        for i, c in enumerate(parent['children']):
            if child_name == c["name"]:
                parent['children'].insert(0, parent['children'].pop(i))
                # del parent['children'][i]
                print(f"Moved {child_name} to bottom")
                return
            else:
                reorder_child(c, child_name)

# %%
# from tree, prune child node and place it as first child of parent node
def move_to_bottom(muts: list[str]):
    for mut in muts:
        res = get_node_with_mutation(auspice["tree"], mut)
        # print(res)
        child=res['child']
        reorder_child(modified["tree"], child['name'])
#%%
# move_to_bottom(['C10507T','C29200T','C9866T'])
move_to_bottom(['C10507T','C9866T'])
# move_to_bottom(['C9866T'])
# res2 = get_node_with_mutation(modified['tree'], 'C10507T')

# parent2=res2['parent']
# child2=res2['child']

# prune_child(modified['tree'], child2['name'])

# res = get_node_with_mutation(auspice['tree'], 'C9866T')
# parent=res['parent']
# child=res['child']
# prune_child(modified['tree'], child['name'])


json.dump(modified, open('auspice/BA.3.2.json', 'w'), indent=2)
# %%

# Now prune

def prune_children(parent, child_name):
    if 'children' in parent:
        for i, c in enumerate(parent['children']):
            if child_name == c["name"]:
                del c["children"]
                c["node_attrs"]["div"] = 11
                print(f"Pruned {child_name} and set div to 11")
                return
            else:
                prune_children(c, child_name)

to_prune = get_node_with_mutation(modified['tree'], 'C897A', max_div=100)
# print(to_prune)
prune_children(modified["tree"], to_prune["child"]["name"])
json.dump(modified, open('auspice/pruned/BA.3.2.json', 'w'), indent=2)



# %%
