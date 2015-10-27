#!/usr/bin/python3.4

"""
Here we set some default values.
They are mostly file names and paths, and default parameters to be given to
the R scripts.
"""

###############################################################################

RCODE     = "/home/rilla/nucleServ/rcode"
RUN_R     = "/usr/bin/Rscript"
LOAD_BAMS = "readBAM.R"
NUCLER    = "nucleR.R"
NUCDYN    = "nucleosomeDynamics.R"
# OUT_DIR   = "/orozco/services/R-data"
# IN_DIR    = "/orozco/services/R-data/in_bams"

###############################################################################

# Default values
DEFAULT_OPTS = {# General options
                "gen":    {"wd":    None,
                           "cores": None,
                           "f1":    None,
                           "f2":    None,
                           "type1": None,
                           "type2": None},

                # Options specific to nucleR
                "nucleR": {"fdrOverAmp":  "0.05",
                           "fragmentLen": None,
                           "components":  "1",
                           "trim":        "50",
                           "pcKeepComp":  "0.02",
                           "width":       "125",
                           "threshold":   "35%",
                           "dyadlength":  None,
                           "minoverlap":  None,
                           "wthresh":     "0.6",
                           "hthresh":     "0.4"},

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
