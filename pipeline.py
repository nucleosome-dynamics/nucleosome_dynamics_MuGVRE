#!/usr/bin/python3.4

"""
Main wrapper script
This will run nucleR and NucDyn on a pair of experiments
"""

###############################################################################

import sys

from helpers import parse_args, get_opts, mkdir_p, get_args_ls
from exp_methods import Experiment, nucleosome_dynamics

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

    actions = {"preproc1": lambda: exp1.load(),
               "preproc2": lambda: exp2.load(),
               "nucleR1":  lambda: exp1.nucleR(cores, *nucleR_optargs),
               "nucleR2":  lambda: exp2.nucleR(cores, *nucleR_optargs),
               "nucdyn":   lambda: nucleosome_dynamics(exp1, exp2,
                                                       wd, cores,
                                                       *nucdyn_optargs)}
    for c in calcs:
        actions[c]()

    return 0

###############################################################################

if __name__ == "__main__":
    config_f, calcs = parse_args()
    sys.exit(main(config_f, calcs))

############################################################################### 
