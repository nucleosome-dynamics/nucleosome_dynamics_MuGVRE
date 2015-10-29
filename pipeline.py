#!/usr/bin/python3.4

"""
Main wrapper script
This will run nucleR and NucDyn on a pair of experiments
"""

###############################################################################

import sys

from helpers import parse_args, get_opts, get_args_ls

from classes.experiment import Experiment
from classes.calculations import Load, NucleR, NucDyn
from classes.run import Run

###############################################################################

def main(config_f, calcs):
    # get the input arguments
    opts = get_opts(config_f)

    cores = opts["gen"]["cores"]

    wd = opts["gen"]["wd"]
    run = Run(wd)

    # wrap the information for each experiment in an Experiment object
    exp1 = Experiment(opts["gen"]["f1"], opts["gen"]["type1"])
    exp2 = Experiment(opts["gen"]["f2"], opts["gen"]["type2"])

    nucleR_optargs = get_args_ls(opts["nucleR"])
    nucdyn_optargs = get_args_ls(opts["NucDyn"])

    load1 = Load(exp1, 1, run)
    load2 = Load(exp2, 2, run)

    nucleR1 = NucleR(exp1, nucleR_optargs, cores, 1, run)
    nucleR2 = NucleR(exp2, nucleR_optargs, cores, 2, run)

    nucDyn = NucDyn(exp1, exp2, nucdyn_optargs, cores, run)

    actions = {"preproc1": load1,
               "preproc2": load2,
               "nucleR1":  nucleR1,
               "nucleR2":  nucleR2,
               "nucdyn":   nucDyn}

    for c in calcs:
        actions[c].go()

    return 0

###############################################################################

if __name__ == "__main__":
    config_f, calcs = parse_args()
    sys.exit(main(config_f, calcs))

############################################################################## 
