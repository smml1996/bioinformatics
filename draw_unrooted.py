import networkx as nx
import matplotlib.pyplot as plt

G=nx.Graph()

file1 = open('../results/tree_matrix.txt', 'r')

lines = file1.readlines()
nodes = int(lines[0]);


chars = 'ABCDEFGHIJKLMNOPQRSTUWXYZ'
c = 0
for i in range(1,len(lines)):
    elements = lines[i].split()
    for j in range(c+1, len(elements)):
        temp = int(elements[j])
        print(temp)
        if temp != -1:
            if c < nodes:
                if j < nodes:
                    G.add_edge(chars[c],chars[j],weight=elements[j])
                else:
                    G.add_edge(chars[c],'N'+ str(j-nodes),weight=elements[j])
            else:
                if j < nodes:
                    G.add_edge('N' + str(c-nodes),chars[j],weight=elements[j])
                else:
                    G.add_edge('N' + str(c-nodes),'N'+ str(j-nodes),weight=elements[j])
    c+=1
pos = nx.spring_layout(G)


nx.draw(G,pos,with_labels = True)
labels = nx.get_edge_attributes(G,'weight')
nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
plt.savefig("../results/graph")