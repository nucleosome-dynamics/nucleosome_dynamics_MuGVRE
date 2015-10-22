#!/usr/bin/python3.4

"""
Miscellanius helper functions
"""

###############################################################################

import os
import subprocess

from defaults import DEFAULT_OPTS

###############################################################################

def mkdir_p(d):
    try:
        os.makedirs(d)
    except FileExistsError:
        pass


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


def parse_exp_name(x):
    return x.split('/')[-1].split('.')[0]


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
