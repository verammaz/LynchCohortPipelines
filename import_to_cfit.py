import argparse
import os
import json

from cfit.util.Analysis import Analysis
from cfit.plot.PlotTree import PlotTree
from cfit.plot.PlotTreeAll import PlotTreeAll

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-hdir', required=True, help="Analysis data home directory")
    parser.add_argument('-patient_id', required=True, help='Patient identifier')
    parser.add_argument('-drivers', nargs='+', help='List of driver mutations to mark on clone tree plot')


    args = parser.parse_args()

    # Read files with driver gene names and recurring antigens

    drivers = args.drivers

    if drivers is None:
        drivers = []

        if os.path.exists(os.path.join(args.hdir, 'adenoma_drivers.txt')):
            with open(os.path.join(args.hdir, 'adenoma_drivers.txt'), 'r') as f:
                drivers_adenoma = [line.strip() for line in f.readlines()]
                drivers.extend(drivers_adenoma)

        if os.path.exists(os.path.join(args.hdir, 'carcinoma_drivers.txt')):
            with open(os.path.join(args.hdir, 'carcinoma_drivers.txt'), 'r') as f:
                drivers_carcinoma = [line.strip() for line in f.readlines()]
                drivers.extend(drivers_carcinoma)

        if os.path.exists(os.path.join(args.hdir, 'adenoma_shared_fs_interpatient.txt')):
            with open(os.path.join(args.hdir, 'adenoma_shared_fs_interpatient.txt'), 'r') as f:
                fs_recur_interpatient = [line.strip() for line in f.readlines()]
                drivers.extend(fs_recur_interpatient)

        if os.path.exists(os.path.join(args.hdir, 'adenoma_shared_fs_intrapatient.txt')):
            with open(os.path.join(args.hdir, 'adenoma_shared_fs_intrapatient.txt'), 'r') as f:
                fs_recur_intrapatient = [line.strip() for line in f.readlines()]
                drivers.extend(fs_recur_intrapatient)


    anl = Analysis()
    anl.set_MHC_version("pan41")
   
    with open(os.path.join(args.hdir, "cfit_config.json")) as f:
        config = json.load(f)

    with open(os.path.join(args.hdir, "cfit_mapping.json")) as f1:
        mappingjs = json.load(f1)


    mapping = []

    # get patient mapping
    for data_dict in mappingjs:
        if data_dict['name'] == args.patient_id or data_dict['pname'] == args.patient_id:
            mapping.append(data_dict)

    # initialize cfit patient(s):
    anl.initialize_config(config, mapping, dir=args.hdir, kd_thr=500, ns=[9], tree_format="pairtree")

    # to generate all graphs and save into single html file:
    for patient in anl.patients.values():
        pat_dir = os.path.join(args.hdir, config["tree_dir"], patient.name)
        plots = PlotTreeAll(patient, drivers=drivers, outdir=pat_dir, tree_format='pairtree')

