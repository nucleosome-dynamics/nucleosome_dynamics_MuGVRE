#!/home/rilla/anaconda/envs/py3k/bin/python

###############################################################################

import os
import subprocess
import sys

###############################################################################

SCRIPTS_DIR = "/home/rilla/scratch/nucleosome_dynamics/webserver/scripts"
BASE_OUT_DIR = "/home/rilla/scratch/nucleosome_dynamics/webserver/tmp"
RUN_R = "/opt/R-2.15.3/bin/Rscript"
LOAD_BAMS = "readBAM.R"
NUCLER = "nucleR.R"
NUCDYN = "nucleosomeDynamics.R"

###############################################################################

# Default values
DEFAULT_OPTS = {# General options
                "gen":    {"cores":          None,
                           "f1":             None,
                           "f2":             None,
                           "type1":          None,
                           "type2":          None},

                # Options specific to nucleR
                "nucleR": {"fdrOverAmp":     "0.05",
                           "fragmentLen":    None,
                           "components":     "1",
                           "trim":           "50",
                           "pcKeepComp":     "0.02",
                           "width":          "125",
                           "threshold":      "35%",
                           "dyadlength":     None,
                           "minoverlap":     None,
                           "wthresh":        "0.6",
                           "hthresh":        "0.4"},

                # Options specific to NucDyn
                "NucDyn": {"maxLen":         "170",
                           "equalSize":      "FALSE",
                           "roundPow":       "5",
                           "readSize":       "140",
                           "maxDiff":        None,
                           "rangeStart":     None,
                           "rangeEnd":       None,
                           "chr":            None,
                           "combined":       "TRUE",
                           "same.magnitude": "2",
                           "threshold":      "60%"}}

###############################################################################

def mkdir_p(d):
    try:
        os.makedirs(d)
    except FileExistsError:
        pass


def read_config(inf, default={"gen":{}, "nucleR":{}, "NucDyn":{}}):
    """
    Read specified options in config_f
    """
    opts = dict(default)
    with open(config_f) as fh:
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

###############################################################################

config_f = sys.argv[1]

out_dir = BASE_OUT_DIR + "/foo"
load_bams, nucleR, nucdyn = ("{0}/{1}".format(SCRIPTS_DIR, i)
                             for i in (LOAD_BAMS, NUCLER, NUCDYN))
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

mkdir_p(out_dir)

exp1, exp2 = (parse_exp_name(opts["gen"][i])
              for i in ("f1", "f2"))

###############################################################################

rs1, rs2 = ("{0}/{1}.RData".format(out_dir, x)
            for x in (exp1, exp2))

#subprocess.call([RUN_R, load_bams,
#                 "--type", opts["gen"]["type1"],
#                 "--input", opts["gen"]["f1"],
#                 "--output", rs1])
#
#subprocess.call([RUN_R, load_bams,
#                 "--type", opts["gen"]["type2"],
#                 "--input", opts["gen"]["f2"],
#                 "--output", rs2])

###############################################################################

nucleR_out1, nucleR_out2 = ("{0}/{1}.gff".format(out_dir, i)
                            for i in (exp1, exp2))
nucleR_optargs = get_args_ls(opts["nucleR"])

subprocess.call([RUN_R, nucleR,
                 "--input", rs1,
                 "--cores", opts["gen"]["cores"],
                 "--output", nucleR_out1,
                 "--type", opts["gen"]["type1"]] + nucleR_optargs)

subprocess.call([RUN_R, nucleR,
                 "--input", rs2,
                 "--cores", opts["gen"]["cores"],
                 "--output", nucleR_out2,
                 "--type", opts["gen"]["type2"]] + nucleR_optargs)

###############################################################################

nucdyn_out = "{0}/{1}_{2}.gff".format(out_dir, exp1, exp2)
nucdyn_optargs = get_args_ls(opts["NucDyn"])

#subprocess.call([RUN_R, nucdyn,
#                 "--input1", rs1,
#                 "--input2", rs2,
#                 "--output", nucdyn_out,
#                 "--cores", opts["gen"]["cores"]] + nucdyn_optargs)

###############################################################################
