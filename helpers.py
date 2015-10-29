#!/usr/bin/python3.4

"""
Miscellanius helper functions
"""

###############################################################################

import argparse
import os
import subprocess
import sys

from defaults import DEFAULT_OPTS

###############################################################################

def parse_args():
    descr = "Wrapper to run nucleR and NucleosomeDynamics"
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument("-c",
                        "--config",
                        metavar="config_file",
                        required=True,
                        help="configuration file")

    parser.add_argument("-a",
                        "--all",
                        action="store_true",
                        help="run the whole pipeline")

    calc_args = parser.add_argument_group('single calculations')

    calc_args.add_argument("--preproc1",
                           action="store_true",
                           help="preprocessing of the first experiment")

    calc_args.add_argument("--preproc2",
                           action="store_true",
                           help="preprocessing of the second experiment")

    calc_args.add_argument("--nucleR1",
                           action="store_true",
                           help="nucleR on the first experiment")

    calc_args.add_argument("--nucleR2",
                           action="store_true",
                           help="nucleR on the second experiment")

    calc_args.add_argument("--nucdyn",
                           action="store_true",
                           help="NucleosomeDynamics")

    args = parser.parse_args()

    avail_calcs = ["preproc1",
                   "preproc2",
                   "nucleR1",
                   "nucleR2",
                   "nucdyn"]

    if args.all:
        todo = avail_calcs
    else:
        todo = [calc
                for calc in avail_calcs
                if getattr(args, calc)]

    if not todo:
        print("nothing to be done")
        sys.exit(0)

    return args.config, todo

###############################################################################



def read_config(inf, default={"gen": {}, "nucleR": {}, "NucDyn": {}}):
    """
    Read specified options in config_f
    """
    opts = dict(default)
    with open(inf) as fh:
        for line in fh:
            vals = line.split()
            if vals[0] in ("nucleR", "NucDyn"):
                group, param, value = vals
            else:
                param, value = vals
                group = "gen"
            opts[group][param] = value
    return opts


def set_by(d, xs, y, f=lambda x: x):
    """
    Some unspeficied options will be set to a value that depends on another
    parameter
    """
    new_dict = dict(d)
    for x in xs:
        if d[x] is None:
            new_dict[x] = f(d[y])
    return new_dict


def rm_nones(d):
    """
    Options that are still None in either nucleR or NucDyn won't be passed
    as arguments
    """
    return {k:v for k, v in d.items() if v is not None}


def check_non_optionals(d):
    """
    All general options have to be specified, so there can be no None
    left there
    """
    unspecified = [k for k, v in d.items() if v is None]
    if unspecified:
        print("Some non-optional parameters were not specified:")
        print("\t-", *unspecified)
        sys.exit(1)


def flatten(xss):
    return [x for xs in xss for x in xs]


def get_args_ls(d):
    return flatten((("--" + k, v) for k, v in d.items()))


def get_opts(config_f):
    """
    Read and check the input optional arguments.
    """
    opts = read_config(config_f, DEFAULT_OPTS)
    opts["nucleR"] = set_by(opts["nucleR"],
                            ["dyadlength", "minoverlap"],
                            "trim")
    opts["NucDyn"] = set_by(opts["NucDyn"],
                            ["maxDiff"],
                            "readSize",
                            lambda x: str(int(x)//2))
    for k in ("NucDyn", "nucleR"):
        opts[k] = rm_nones(opts[k])
    check_non_optionals(opts["gen"])
    return opts

###############################################################################
