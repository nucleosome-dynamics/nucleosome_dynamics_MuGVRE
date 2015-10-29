#!/usr/bin/python3.4

"""
Main wrapper script
This will run nucleR and NucDyn on a pair of experiments
"""

###############################################################################

import sys

from helpers import parse_args, get_opts, mkdir_p, get_args_ls
from exp_methods import Experiment, Load, NucleR, NucDyn

###############################################################################

def main(config_f, calcs):
    # get the input arguments
    opts = get_opts(config_f)

    wd = opts["gen"]["wd"]
    cores = opts["gen"]["cores"]

    # create the working directory where files will be stored
    mkdir_p(wd)

    # wrap the information for each experiment in an Experiment object
    exp1 = Experiment(opts["gen"]["f1"], opts["gen"]["type1"], wd)
    exp2 = Experiment(opts["gen"]["f2"], opts["gen"]["type2"], wd)

    nucleR_optargs = get_args_ls(opts["nucleR"])
    nucdyn_optargs = get_args_ls(opts["NucDyn"])

    info_file = wd + "/info.txt"

    load1 = Load(exp1, 1, info_file)
    load2 = Load(exp2, 2, info_file)

    nucleR1 = NucleR(exp1, nucleR_optargs, cores, 1, info_file)
    nucleR2 = NucleR(exp2, nucleR_optargs, cores, 2, info_file)

    nucDyn = NucDyn(exp1, exp2, wd, nucdyn_optargs, cores, info_file)

    actions = {"preproc1": load1,
               "preproc2": load2,
               "nucleR1":  nucleR1,
               "nucleR2":  nucleR2,
               "nucdyn":   nucDyn}

    for c in calcs:
        actions[c].run()

    return 0

###############################################################################

if __name__ == "__main__":
    config_f, calcs = parse_args()
    sys.exit(main(config_f, calcs))

############################################################################## 
