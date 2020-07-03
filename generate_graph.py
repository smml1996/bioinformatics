import networkx as nx
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import matplotlib.pyplot as plt

G=nx.Graph()

file1 = open('../results/graph.txt', 'r')

lines = file1.readlines()
nodes = int(lines[0]);
for i in range(len(lines) -1, 0, -1):

    elements = lines[i].split()
    G.add_edge(elements[0],elements[1],weight=elements[2])

pos=graphviz_layout(G, prog='dot')
nx.draw(G,pos,with_labels = True)
labels = nx.get_edge_attributes(G,'weight')
nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
plt.savefig("../results/graph")