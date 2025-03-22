import argparse
import os
import json
from collections import defaultdict

from cfit.util.Analysis import Analysis
from cfit.plot.PlotTree import PlotTree


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-hdir', required=True, help="Analysis data home directory")
    parser.add_argument('-patient_id', required=True, help='Patient/Sample identifier')
    parser.add_argument('-gene_ann_files', nargs='+', help='File(s) with list of driver mutations to mark on clone tree plot')
    parser.add_argument('-config', required=False, default='cfit_config.json')
    parser.add_argument('-mapping', required=False, default='cfit_mapping.json')
    parser.add_argument('-outdir', required=False)

    args = parser.parse_args()

    # Read files with driver gene names and recurring antigens

    gene_ann_files = args.gene_ann_files

    drivers = []

    if gene_ann_files is not None:
        for file in gene_ann_files:
            with open(os.path.join(args.hdir, file), 'r') as f:
                for line in f.readlines():
                    drivers.append(line.strip())

    
    anl = Analysis()

    with open(os.path.join(args.hdir, args.config)) as f1:
        config = json.load(f1)

    with open(os.path.join(args.hdir, args.mapping)) as f2:
        mappingjs = json.load(f2)


    mapping = []

    # get patient mapping:
    for data_dict in mappingjs:
        if data_dict['name'] == args.patient_id or data_dict['pname'] == args.patient_id:
            mapping.append(data_dict)

    assert(len(mapping) == 1)

    # initialize cfit:
    anl.initialize_config(config, mapping, dir=args.hdir, kd_thr=500, ns=[9], tree_format="pairtree")

    assert(len(anl.patients.values()) == 1)
    patient = set(anl.patients.values()).pop()

    outdir = os.path.join(args.hdir, config["tree_dir"], patient.name)
    
    plottree = PlotTree(patient, 'topology', drivers, 0)
    plottree.plot()
    
    filename = f"{patient.name}_clonetree.html"

    with open(os.path.join(outdir, filename), 'w') as f:
        f.write("<html> <body> <h1> Clone Tree for \
        <font color = #000000>{}</font></h1>\n \
        </body></html>".format(patient.name))
        f.write(plottree.figs[0].to_html(full_html=False, include_plotlyjs='cdn'))
    

