import networkx as nx
from functools import reduce
import mpbn

### From boolean network to influence graph
def influence_graph(boolean_model):
    if boolean_model.split(".")[1] == "dot":
        return(nx.nx_agraph.read_dot(boolean_model))
    if boolean_model.split(".")[1] == "bnet":
        f = mpbn.load(boolean_model) #./min-1.bnet
        return(f.influence_graph())
#todo : add the case .bn

### Misc 
def input_nodes(G):
    input_nodes=[]
    for node in list(G.nodes):
        if list(G.predecessors(node)) == [] and not(list(G.successors(node)) == []):
                input_nodes.append(node)
    if input_nodes == []:
        print("No input nodes")
    return input_nodes

def output_nodes(G):
    output_nodes=[]
    for node in list(G.nodes):
        if list(G.successors(node)) == [] and not(list(G.predecessors(node)) == []):
                output_nodes.append(node)
    if output_nodes == []:
        print("No output nodes")
    return output_nodes

def add_path(G, subset_nodes, list_paths):
    G_sub=G.subgraph(subset_nodes).copy()
    for path in list_paths:
        #print(path[::len(path)-1])
        nx.add_path(G_sub, path[::len(path)-1], color='red')
        return(G_sub)

### About stongly connected components
def info_strongly_connected_components(G):
    sccs = [list(scc) for scc in nx.strongly_connected_components(G) if len(scc) == 1 and list(scc)[0] in list(nx.nodes_with_selfloops(G))]
    print(f"Number of self-looped nodes that are a strongly connected components : {len(sccs)}")
    sccs = [list(scc) for scc in nx.strongly_connected_components(G) if len(scc) > 1]
    print(f"Number of strongly connected components (no singleton) : {len(sccs)}")
    return(sccs)

def subset_node_sccs(G):
    feedback = reduce(set.union, 
                  (set(scc) 
                   for scc in nx.strongly_connected_components(G) 
                   if len(scc) > 1))
    return(G.subgraph(feedback), feedback)
#todo : enlever le 0, sert à rien. 

### About supplementary biomarkers
def intersection(G,biomarkers):
    marker_in_G=[]
    for marker in biomarkers:
        if G.has_node(marker):
            marker_in_G.append(marker)
    return(marker_in_G)

### About the paths between SCCs and downstream biomarkers
#todo upstream paths, only downstream markers so far 

def list_markers_not_in_G(G, markers):
    for n1 in markers:
        if n1 not in G: 
            print(f"Marker {n1} not in the influence graph")

def info_path_markers_sccs(G,sccs,markers,direction):
    upstream_paths={} #marker gene: (scc,scc_gene)
    downstream_paths={}
    for scc in sccs:
        for n1 in markers:
            for n2 in scc:
                if n1 in G: 
                    if direction=="upstream" and nx.has_path(G,n1,n2):
                        upstream_paths[n1]=(str(scc),n2)
                    if direction=="downstream" and nx.has_path(G,n2,n1):
                        downstream_paths[(str(scc),n2)]=(n1)
    if direction == "upstream":
        if len(upstream_paths)==0:
            print("No upstream markers")
            return False
        else:
            print("Upstream markers: checked")
            return upstream_paths
    if direction=="downstream":
        if len(downstream_paths)==0:
            print("No downstream markers")
            return False
        else:
            print("Downstream markers: checked")
            return downstream_paths
        
def info_markers(markers, G, sccs):
    info_path_markers_sccs(G,sccs,markers,"downstream")
    info_path_markers_sccs(G,sccs,markers,"upstream")

        
def shortest_path_marker_scc(G, marker, scc):
    feedback = subset_node_sccs(G) #feedback[1] is the set of nodes within a scc
    shortest_path_among_scc=[]
    last_len_shortest_path=10000000 #high to begin with, need a proper initialization (but flemme)
    for node_scc in scc:
        ### if marker is in a strongly connected component, no path to search, already treated. 
        if marker not in feedback[1] and marker in G:
            ### downtsream markers
            if nx.has_path(G,node_scc,marker): #(here sorting of downstream)
                #print("downstream marker :", marker)
                shortest_path_scc=nx.shortest_path(G, source=node_scc, target=marker)
                if (len(shortest_path_scc)<last_len_shortest_path):
                    last_len_shortest_path=len(shortest_path_scc)
                    shortest_path_among_scc=shortest_path_scc
            ### upstream markers
            if nx.has_path(G,marker, node_scc) : #(here sorting of upstream)
                #print("upstream marker :", marker)
                shortest_path_scc=nx.shortest_path(G, source=marker, target=node_scc)
                if (len(shortest_path_scc)<last_len_shortest_path):
                    last_len_shortest_path=len(shortest_path_scc)
                    shortest_path_among_scc=shortest_path_scc
    return(shortest_path_among_scc)

def list_paths(G, markers_in_G, sccs): #list_path
    downstream_path_sccs=[] #path
    for marker in markers_in_G:
        for scc in sccs:
            downstream_path_sccs.append(shortest_path_marker_scc(G, marker,scc)) 
    return(list(filter(None, downstream_path_sccs)))  

### About the paths between the SCCs
def shortest_path_sccs(G,H):
    path_between_sccs=[]
    if list(H.edges()) == []:
        print(f"No path between the {len(H)} strongly connected components")
    for edge in H.edges():
        last_len_shortest_path=100000
        shortest_path=[]
        scc_source_edges=list(H.nodes.data()[int(edge[0])]['members'])
        scc_target_edges=list(H.nodes.data()[int(edge[1])]['members'])
        for node1 in scc_source_edges:
            for node2 in scc_target_edges:
                path=nx.shortest_path(G, source=node1, target=node2)
                if len(path)<last_len_shortest_path:
                    last_len_shortest_path=len(path)
                    shortest_path=path
        path_between_sccs.append(shortest_path)
    return(path_between_sccs)

def draw_compressed_graph(G, markers, sccs, outputfile):
    feedback = subset_node_sccs(G)
    paths = list_paths(G, markers, sccs)
    G_compressed=G.subgraph(feedback[1]).copy()
    ### path between sccs and markers
    for path in paths:
        nx.add_path(G_compressed, path[::len(path)-1], color='red')
    ### path between sccs
    H = nx.condensation(G.subgraph(feedback[1]).copy())
    path_between_sccs=shortest_path_sccs(G,H)
    G_compressed.remove_edges_from(path_between_sccs)
    for path in path_between_sccs:
        nx.add_path(G_compressed, path, color='blue')
    ### Save in .dot
    try :
        nx.drawing.nx_pydot.write_dot(G_compressed,outputfile)
        print("The compression of the influence graph is done.")
        print(f"Please run the following command in your terminal : dot -Tpdf {outputfile}.dot -o {outputfile}.pdf")
        print("#### END ####")
    except FileNotFoundError : 
        print("The path of the file is not correct.")