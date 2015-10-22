#!/usr/bin/python3.4

"""
Main wrapper script
This will run nucleR and NucDyn on a pair of experiments
"""

###############################################################################

import sys

from defaults import OUT_DIR
from helpers import get_opts, mkdir_p, get_args_ls
from exp_methods import Experiment, nucleosome_dynamics

############################################################################### 

def main(config_f):

    # get the input arguments
    opts = get_opts(config_f)

    # create the working directory where files will be stored
    wd = "{0}/{1}".format(OUT_DIR, opts["gen"]["wd"])
    mkdir_p(wd)

    # wrap the information for each experiment in an Experiment object
    exp1, exp2 = (Experiment(opts["gen"]["f{0}".format(i)],
                             opts["gen"]["type{0}".format(i)],
                             wd)
                  for i in (1, 2))

    # run nucleR
    nucleR_optargs = get_args_ls(opts["nucleR"])
    for e in (exp1, exp2):
        e.nucleR(opts["gen"]["cores"], *nucleR_optargs)

    # run NucDyn
    nucdyn_optargs = get_args_ls(opts["NucDyn"])
    nucleosome_dynamics(exp1, exp2, wd, opts["gen"]["cores"], *nucdyn_optargs)

    return 0

###############################################################################

if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1]))
    except IndexError:
        print("An input configuration file must be supplied")
        sys.exit(1)

###############################################################################
