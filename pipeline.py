usr/bin/python3.4

###############################################################################

import subprocess
import sys
from defaults import *
from helpers import *

###############################################################################

config_f = sys.argv[1]

load_bams, nucleR, nucdyn = ("{0}/{1}".format(RCODE, i)
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

mkdir_p(OUT_DIR)
mkdir_p(tmp_dir)

exp1, exp2 = (parse_exp_name(opts["gen"][i])
              for i in ("f1", "f2"))

###############################################################################

rs1, rs2 = ("{0}/{1}.RData".format(tmp_dir, x)
            for x in (exp1, exp2))

inf1, inf2 = ("{0}/{1}".format(IN_DIR, opts["gen"][x])
              for x in ("f1", "f2"))

subprocess.call([RUN_R, load_bams,
                 "--type", opts["gen"]["type1"],
                 "--input", inf1,
                 "--output", rs1])

subprocess.call([RUN_R, load_bams,
                 "--type", opts["gen"]["type2"],
                 "--input", inf2,
                 "--output", rs2])

###############################################################################

nucleR_out1, nucleR_out2 = ("{0}/{1}.gff".format(OUT_DIR, i)
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

nucdyn_out = "{0}/{1}_{2}.gff".format(OUT_DIR, exp1, exp2)
nucdyn_optargs = get_args_ls(opts["NucDyn"])

subprocess.call([RUN_R, nucdyn,
                 "--input1", rs1,
                 "--input2", rs2,
                 "--output", nucdyn_out,
                 "--cores", opts["gen"]["cores"]] + nucdyn_optargs)

###############################################################################
