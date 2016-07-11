import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy
import igraph
import itertools
import pickle


#Always inserts new nodes to the right of the focal node
#Split will create a new identical node and insert to the right of the focal node
class SplitListNode:
    def __init__(self, ref_object):
        self.ref_object = ref_object
        self.left_node = None
        self.right_node = None
        self.split_clones = []

    def split_node(self):
        temp_left = self.left_node
        temp_right = self.right_node

        new_node = SplitListNode(self.ref_object)
        new_node.left_node = self
        new_node.right_node = temp_right
        self.right_node = new_node

        new_node.split_clones.append(self)
        self.split_clones.append(new_node)

        return new_node

    def insert_node(self, node_to_insert):
        temp_right = self.right_node
        node_to_insert.left_node = self
        node_to_insert.right_node = temp_right
        self.right_node = node_to_insert

class SplitList:
    def __init__(self, root_object=None):
        self.nodes = [SplitListNode(root_object)]

    def split_node(self, node):
        new_node = node.split_node()
        self.nodes.append(new_node)
        return new_node

    def split_and_insert(self, node_to_split, node_to_insert):
        #TODO
        #Actually want to look at all split nodes belonging to same class, and decide which
        #to split again -- to lay ou mutations better
        split_node = self.split_node(node_to_split)
        self.nodes.append(node_to_insert)
        node_to_split.insert_node(node_to_insert)
        return split_node

    def __str__ (self):
        return str([str(j) for j in self.nodes])

    def get_start_node(self):
        for n in self.nodes:
            if n.left_node == None:
                return n



class MullerFisherList:
    def __init__(self, num_generations):
        self.node_counter = itertools.count(0)
        self.graph = igraph.Graph(directed = True,
                                 graph_attrs = {'population_size': 1,
                                                'generations': num_generations},
                                 vertex_attrs = {'frequencies': np.ones(num_generations),
                                                'split_list_node': None,
                                                'plot_color': 'r'},
                                 edge_attrs = {})

        self.graph.add_vertex(name=next(self.node_counter), frequencies=np.ones(num_generations), plot_color=np.random.rand(3,1))
        
        #how meta!
        self.split_list = SplitList(self.graph.vs[0])
        self.graph.vs[0]['split_list_node'] = self.split_list.nodes[0]

    #only have to update parent's freq!
    def add_mutation(self, parent_mut_node, new_mut_frequencies):
        self.graph.add_vertex(name=next(self.node_counter), frequencies=new_mut_frequencies, plot_color=np.random.rand(3,1))
        self.graph.add_edge(parent_mut_node, self.graph.vs[self.graph.vcount()-1])

        new_splitlist_node = SplitListNode(self.graph.vs[self.graph.vcount()-1])
        self.graph.vs[self.graph.vcount()-1]['split_list_node'] = new_splitlist_node

        split_node = self.split_list.split_and_insert(parent_mut_node['split_list_node'], new_splitlist_node)
        #handle updating of frequency data here.
        #get split clones, divide freq vector by total number of splits (is that it?!)
        num_splits = len(split_node.split_clones) + 1
        split_node.ref_object['frequencies'] = split_node.ref_object['frequencies'] - new_mut_frequencies

        #divide freq by number of splits we have
        split_node.ref_object['frequencies'] = split_node.ref_object['frequencies'] / float(num_splits)



NUM_CYCLES = 300

def make_freq_vect(node, num_cycles=NUM_CYCLES):
    yvec = [0 for i in range(0, num_cycles)]
    for j in range(len(node['abundances'])):
        yvec[j+node['first_seen']] = node['abundances'][j]
    return yvec



#genotypes = pickle.load(open('tree-end.pickle', 'rb'))
mutation_nodes = genotypes.vs.select(genotype_node_eq = False)
genotype_nodes = genotypes.vs.select(genotype_node_eq = True)

a = MullerFisherList(NUM_CYCLES)
#print(genotypes.neighbors(genotypes.vs[0], mode="out"))
#dfs_tuple = genotypes.subgraph(genotype_nodes).bfs(0,mode="out")
print(list(genotypes.vs[genotypes.neighbors(genotype_nodes[1], mode="in")]))




'''

for mut_node in mutation_nodes:
    plt.plot(make_freq_vector(mut_node['abundances'], mut_node['first_seen'], mut_node['last_seen']))

plt.show()
print(mutation_nodes[0])

print(genotypes.neighbors(mutation_nodes[0]))


a = MullerFisherList(100)
a.add_mutation(a.graph.vs[0], y_1)
a.add_mutation(a.graph.vs[1], y_2/4)

start_node = a.split_list.get_start_node()
to_plot = []
to_plot_colors = []

next_node = start_node
while next_node != None:
    to_plot.append(next_node.ref_object['frequencies'])
    to_plot_colors.append(next_node.ref_object['plot_color'])
    next_node = next_node.right_node


fig = plt.figure()
ax1 = fig.add_subplot(111)
'''
