import pandas as pd
import sys
sys.path.append('../')
import graph_analysis as ga

#### About pĥenotypic markers ####
#############################################

#From PNAS 2023 Figure 3
ST_HSC = ["Pdzk1ip1", "Apoe", "Trim47", "ADT_CD150"]
Progenitors =  ["Adgrl4", "Car2", "Kit", "Mpo", "Ctsg", "Igfbp4", "Bcl2", "Dach1", "Lyz2", "Ccr2", "ADT_CD117", "ADT_CD115"]
Lymphoid_genes = ["Dntt", "Rag2", "Myl10", "Ebf1"]
Dendritic_cell_development = ["Il7r", "ADT_127", "Ly6d", "Irf8", "Cx3cr1", "Bcl11a", "Runx2", "Cd300c", "Clec9a", "Irf4", "Siglech", "Tcf4", "Batf3", "Zeb2", "Id2", "Zbtb46", "ADT_CD105", "ADT_CD11b", "ADT_CD11c", "ADT_CD62L"] 
Cell_cycle = ["Mki67", "Top2a"]
#in gene name, always _ as separator not space ! 

#From Sarah 
ADT_gene = pd.read_csv('./ADT.csv', names=["ADT", "gene"])  
ADTs = [ADT.strip().replace(" ", "_") for ADT in list(ADT_gene["ADT"])]
ADT_gene["ADT"] = ADTs
ADT_gene = ADT_gene.set_index('ADT')
#test : print(ADT_gene.loc['ADT_CD41', 'gene']) #get via index 

marker_category = [ST_HSC, Progenitors, Lymphoid_genes, Dendritic_cell_development, Cell_cycle]
for marker in marker_category:
    for g in marker: 
        try:
            if g[:3] == "ADT":
                #test : print(g, ADT_gene.loc[g, 'gene'].split(", "))
                marker.remove(g)
                marker.append(ADT_gene.loc[g, 'gene'])
        except KeyError:
            print(f"{g} is unknown by Cathy and Sarah.")

#### EXAMPLE HEMAMOX MIN MODEL VERSION 1 ####
#############################################

### Version of the model (here min v1)
# FLAG : I don't know why a relative path doesn't work above... to be investigated. 
hemamox_model = "/home/dboyenval/Documents/bonesis-tools/grntools/grn/example/hemamox_model_version1/min-1.dot" 
G_hemamox = ga.influence_graph(hemamox_model)
print(G_hemamox)

### About markers intersection with the influence graph

hemamox_phenotypic_markers=list(set(ga.intersection(G_hemamox,ST_HSC)
     +ga.intersection(G_hemamox,Progenitors)
     +ga.intersection(G_hemamox,Lymphoid_genes)
     +ga.intersection(G_hemamox,Dendritic_cell_development)
     +ga.intersection(G_hemamox,Cell_cycle)))

print(f"List of markers: {hemamox_phenotypic_markers}")

### About SCCs
sccs_hemamox = ga.info_strongly_connected_components(G_hemamox)

### Upstream markers ? Downstream markers ? ... with regard to SCCs
ga.info_markers(hemamox_phenotypic_markers, G_hemamox, sccs_hemamox)

### Paths between SCCs (to check): 
#H = nx.condensation(G_hemamox.subgraph(feedback[1]).copy())
#path_between_sccs=ga.shortest_path_sccs(G_hemamox,H)
#print(path_between_sccs)

### Drawing the compressed graph 
output = "/home/dboyenval/Documents/bonesis-tools/grntools/grn/example/hemamox_model_version1/G_compressed_hemamox_version1"
ga.draw_compressed_graph(G_hemamox, hemamox_phenotypic_markers, sccs_hemamox, output)

#### EXAMPLE HEMAMOX MIN MODEL VERSION 2 ####
#############################################

### Version of the model (here min v2)
hemamox_model = "/home/dboyenval/Documents/bonesis-tools/grntools/grn/example/hemamox_model_version2/min-1.dot"
G_hemamox = ga.influence_graph(hemamox_model)

sccs_hemamox = ga.info_strongly_connected_components(G_hemamox)
ga.info_markers(hemamox_phenotypic_markers, G_hemamox, sccs_hemamox)

output = "/home/dboyenval/Documents/bonesis-tools/grntools/grn/example/hemamox_model_version2/G_compressed_hemamox_version2"
ga.draw_compressed_graph(G_hemamox, hemamox_phenotypic_markers, sccs_hemamox, output)