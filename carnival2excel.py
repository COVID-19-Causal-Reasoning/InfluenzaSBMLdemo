import argparse
from pathlib import Path
import pandas as pd
import pygraphviz as pg

def export_network(g, activity, nodetype, edgetype, omnipath):
    node_df = pd.DataFrame(dict(activity=[activity[node.attr['fillcolor']]
                                          for node in g],
                                nodetype=[nodetype[node.attr['shape']]
                                          if 'shape' in node.attr 
                                          else 'InferredNode'
                                          for node in g]),        
                            index=g.nodes())
    edge_df = pd.DataFrame([(edge[0],edgetype[edge.attr['arrowhead']], edge[1])
                           for edge in g.edges()],
                          columns=['Source', 'Relationship', 'Target'])
    edge_omnipath = pd.merge( edge_df, omnipath, 
                              how='left',  
                              left_on=['Source','Target'], 
                              right_on=['source_genesymbol','target_genesymbol'] )
    
    return node_df, edge_omnipath


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument( "--perturbation_file", help="Perturbation File", default="NoInput")
    parser.add_argument( "--variable", help="Intervention (e.g. VN1203)")
    parser.add_argument( "--baseline", help="The thing to be compared against. E.g. Mock")
    parser.add_argument( "--constant", help="The things that are held constant for both variable and baseline. For example timepoint (7h) or strain (VN1203")
    parser.add_argument( "--inputdir", help="input directory")
    parser.add_argument( "--outputdir", help="output directory")
    parser.add_argument("--omnipathdb",help="Omnipath file")
    argv = parser.parse_args(["--inputdir", "ResultsCARNIVAL_VN1203vsMock_for_24hr",
                   "--variable", "VN1203", 
                    "--baseline", "Mock",
                   "--constant", "24h",
                    "--omnipathdb", "OmnipathSignedDirectedInteractions.csv",
                   "--outputdir", "ResultsCARNIVAL_VN1203vsMock_for_24hr"])


    #expdir = Path(f'ResultsCARNIVAL_{argv.variable}vs{argv.baseline}_for_{argv.constant}r')
    inputdir = Path(argv.inputdir)
    outputdir = Path(argv.outputdir)
    gv = pg.AGraph(str(expdir/'network_solution.dot'))
    activity = {'mistyrose':'DownRegulated', 'lavender':'UpRegulated'}
    nodetype = {'doublecircle':'FootprintObserved',
                'invhouse':'Perturbation',
                '':'InferredNode'}
    edgetype = {"tee":"inhibits",
                "vee": "activates"}
    omnipath = pd.read_csv(argv.omnipathdb).drop_duplicates()
    nodes, edges = export_network(gv, activity, nodetype, edgetype, omnipath)
    outfile = outputdir/f'{argv.variable}vs{argv.baseline}_for_{argv.constant}r_{argv.perturbation_file}.xlsx' 
    
    with pd.ExcelWriter(str(outfile)) as writer:  
        nodes.to_excel(writer, sheet_name='Protein Activity')
        edges.to_excel(writer, sheet_name='RegulatoryInteractions')
