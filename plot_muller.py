import numpy as np
from matplotlib import pyplot as plt
from copy import deepcopy

# just create some random data
fnx = lambda : np.random.randint(5, 10, 100)
y = np.row_stack((fnx(), fnx(), fnx()))   
# this call to 'cumsum' (cumulative sum), passing in your y data, 
# is necessary to avoid having to manually order the datasets
x = np.arange(100) 

y_0 = [10 for i in range(100)]
y_1 = [np.exp(i/70.0)-1 for i in range(100)]
y_2 = np.zeros(100)
for i in range(80,100):
	y_2[i] = np.sin((i-80)/8.0)



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

	def __str__ (self):
		return str( {"ref": self.ref_object, "left": self.left_node, "right": self.right_node} )

class SplitList:
	def __init__(self, root_object=None):
		self.nodes = [SplitListNode(root_object)]

	def split_node(self, node):
		new_node = node.split_node()
		self.nodes.append(new_node)
		return new_node

	def split_and_insert(self, node_to_split, node_to_insert):
		split_node(node_to_split)
		self.nodes.append(node_to_insert)
		node_to_split.insert_node(node_to_insert)


	def __str__ (self):
		return str([str(j) for j in self.nodes])


a = SplitList(root_object={'frequency':np.ones(100)})

print(a)

'''


def add_decendent_dynamics(parent, offspring):
	parent_half = np.asarray(parent)/2.0
	offspring_sum = np.asarray(np.sum(np.asmatrix(offspring), axis=0))[0]
	#print (offspring_sum)
	parents = np.asarray([parent_half[i] - offspring_sum[i]/2.0 for i in range(len(parent_half))])
	return [parents, offspring, parents]




print(to_plot)

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.stackplot(x, to_plot, baseline="sym")
plt.show()

y_stack = np.cumsum([y_0, y_2, y_1], axis=0)   # a 3x10 array


fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.stackplot(x, y_0, y_2, y_1, colors=['#377EB8','#55BA87','#7E1137'], baseline="sym")
plt.show()
'''
