from rdkit import Chem
import os
import re

import networkx as nx

import matplotlib.pyplot as plt
# parameters part modified base on case

#re_PG_part = re.compile(r'_PG_\(oh~p!*(\%?\d*),oh~s!*(\%?\d*)\)')
#re_TMP_part = re.compile(r'_TMP_\(oh~p!*(\%?\d*),oh~p!*(\%?\d*),oh~p!*(\%?\d*)\)')
#re_ISOPS_part = re.compile(r'_ISOPS_\(cooh!*(\%?\d*),cooh!*(\%?\d*)\)')

re_PG_part = re.compile(r'PG\(oh!*(\%?\d*),oh!*(\%?\d*)\)')
re_TMP_part = re.compile(r'TMP\(oh!*(\%?\d*),oh!*(\%?\d*),oh!*(\%?\d*)\)')
re_ISOPS_part = re.compile(r'ISOPS\(cooh!*(\%?\d*),cooh!*(\%?\d*)\)')



_ISOPS_smile='O=C(O[*])c1cccc(C(=O)O[*])c1'
_PG_smile = 'CC([*])C[*]'
_TMP_smile = 'CCC(C[*])(C[*])C[*]'

_ISOPS_func_index0=[7,24]
_TMP_func_index0=[7,13,18]
_PG_func_index0=[5,10]

def insert_label_cooh(s, index, label):
    if index <0 or index > len(s):
        return s
    if label=='' or label==None:
        # remove [*]
        return s[:index-2] + s[index+1:]
        #return s
    return s[:index] + ':'+str(label) + s[index:]


def insert_label_oh(s, index, label):
    if index <0 or index > len(s):
        return s
    if label=='' or label==None:
        # remove [*]
        return s[:index-2] + 'O'+ s[index+1:]
        #return s
    return s[:index] + ':'+str(label) + s[index:]

def gen_isops_smile_part(isops, fun_index0, bngl_index):
    # insert from end to begin
    s = isops
    if(fun_index0[0] < fun_index0[1]):
        for index, label in reversed(list(zip(fun_index0, bngl_index))):
            s = insert_label_cooh(s, index, label)
    return s

def gen_tmp_smile_part(isops, fun_index0, bngl_index):
    # insert from end to begin
    s = isops
    if(fun_index0[0] < fun_index0[1]):
        for index, label in reversed(list(zip(fun_index0, bngl_index))):
            s = insert_label_oh(s, index, label)
    return s

def gen_pg_smile_part(isops, fun_index0, bngl_index):
    # insert from end to begin
    s = isops
    if(fun_index0[0] < fun_index0[1]):
        for index, label in reversed(list(zip(fun_index0, bngl_index))):
            s = insert_label_oh(s, index, label)
    return s

def bngl_convert_smiles(bngl_frags):
    ss=None
    frags = bngl_frags.split('.')
    linkers = []
    node_attribs=[]
    color_attr = []
    for idx in range(len(frags)):
        unit_frag = frags[idx]

        if 'ISOPS' in unit_frag:
            res=re.findall(re_ISOPS_part, unit_frag)
            fun_index = _ISOPS_func_index0
            smile_temp = _ISOPS_smile
            new_frag = gen_isops_smile_part(smile_temp, fun_index, res[0])
            color_attr.append('blue')
            node_attribs.append('ISOPS')
        elif 'TMP' in unit_frag:
            res = re.findall(re_TMP_part, unit_frag)
            fun_index = _TMP_func_index0
            smile_temp = _TMP_smile
            new_frag = gen_tmp_smile_part(smile_temp, fun_index, res[0])
            color_attr.append('red')
            node_attribs.append('TMP')
        elif 'PG' in unit_frag:
            res = re.findall(re_PG_part, unit_frag)
            fun_index = _PG_func_index0
            smile_temp = _PG_smile
            new_frag = gen_pg_smile_part(smile_temp, fun_index, res[0])
            color_attr.append('yellow')
            node_attribs.append('PG')
        else:
            print("not right bnglfrags")
        #new_frag = gen_smile_part(smile_temp, fun_index, res[0])

        linkers.append(res[0])
        if idx ==0:
            ss = new_frag
        else:
            ss = ss + '.'+ new_frag
    return ss, linkers, node_attribs, color_attr

def find_empty_indices(lst, empty_values=(None, '', [], {}, set())):
    return [index for index, element in enumerate(lst) if element in empty_values]

def bngl_convert_graph(linkers, node_attribs):
    total_node = len(linkers)
    G= nx.Graph()
    active_nodes=[]
    for index, link in enumerate(linkers):
        node = index+1
        unpaired = find_empty_indices(link, empty_values=(''))
        if(len(unpaired)>0):
            active_nodes.append([index, unpaired, node_attribs[index]])
        linked_nodes = [n for n in link if n]
        if(linked_nodes):
            G.add_node(node, name=node_attribs[index])
            for linked_node in linked_nodes:
                for in2, lk2 in enumerate(linkers):
                    if linked_node in lk2 and in2 != index:
                        G.add_edge(node, in2+1)
                        break
    return G, active_nodes



def get_polysmile(bngl_ss):
    mol = Chem.MolFromSmiles(bngl_ss)
    mol_com = Chem.molzip(mol)
    newsmi = Chem.MolToSmiles(mol_com)
    #newsmi = Chem.CanonSmiles(newsmi)
    return newsmi


all_species=None
#with open('./example1_reference.species') as f:
with open('./input_10.species') as f:
    l = f.readlines()
    for i in range(len(l)):
        l[i] = l[i].split()
        if l[i][0][0] == '#':
            l[i] = []
    for i in range(len(l)-1,-1,-1):
        if l[i] == []:
            l.pop(i)
    all_species = l

    for i in range(len(l)):
        specie=l[i][0]
        ss,  edges, node_attribs, color_attr = bngl_convert_smiles(specie)
        print(edges)
        print(specie)
        G, unpaired_node=bngl_convert_graph(edges, node_attribs)
        try: 
            cycle = nx.find_cycle(G)
        except:
            cycle = []
        print(cycle)
        nx.draw(G, node_color=color_attr, with_labels=True)
        plt.show()
        #name = nx.get_node_attributes(G, 'color')
        #print(name)
        smi = get_polysmile(ss)
        print(smi)

#for specie in all_species:
    #print(specie[0])
#    res=re.findall(re_ISOPS_part, specie[0])
    #print(res)

#def change_string(s, pattern, labels):
#    counter =0
#    def replace(match):
#        nonlocal counter
#        label = labels[counter] if counter < len(labels) else ''
#        counter +=1
#        return pattern + label
#    new_s = re.sub(pattern, replace, s)
#    return new_s


#specie='_ISOPS_(cooh!1,cooh)._TMP_(oh~p!2,oh~p!1,oh~p!3)._ISOPS_(cooh,cooh!2)._ISOPS_(cooh,cooh!3)'
#specie='_ISOPS_(cooh!1,cooh)._PG_(oh~p!1,oh~s!2)._ISOPS_(cooh!2,cooh)'
#res = re.findall(re_TMP_part, specie)
#print(res)

#frags = specie.split('.')
#print("BNGL-species: ", frags)
