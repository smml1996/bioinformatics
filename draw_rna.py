import networkx as nx
import matplotlib.pyplot as plt

G=nx.Graph()

file1 = open('../results/complements.txt', 'r')

lines = file1.readlines()
sequence = lines[0]

edge_color = []

for i in range(0, len(sequence)-2):
    G.add_edge(sequence[i]+ str(i),sequence[i+1]+ str(i+1), edge_color='b')
    edge_color.append('b')


for i in range(1, len(lines)):
    elements = lines[i].split()
    edge_color.append('r')
    G.add_edge(sequence[int(elements[0])] + elements[0],sequence[int(elements[1])] + elements[1], edge_color='r')


pos = nx.kamada_kawai_layout(G)
edge_color= nx.get_edge_attributes(G,'edge_color').values()
nx.draw(G,pos,with_labels = True,node_size=500, edge_color=edge_color)
plt.savefig("../results/rna")