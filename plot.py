#!/usr/bin/python3.4

"""
Wrapper to the plotting function
"""

import argparse
import sys

from helpers import get_opts

from classes.experiment import Experiment
from classes.calculations import Plot
from classes.run import Run

def parse_args():
    descr = "Wrapper to the R script for plotting"
    parser = argparse.ArgumentParser(description=descr)

    parser.add_argument("-c",
                        "--config",
                        metavar="config_file",
                        required=True,
                        help="configuration file")

    parser.add_argument("--chr",
                        metavar="chr",
                        required=True,
                        help="chromosome to plot")

    parser.add_argument("--range",
                        metavar="range",
                        required=True,
                        help="range to plot in the chromosome")

    args = parser.parse_args()
    start, end = map(int, args.range.split("-"))

    return args.config, start, end, args.chr


def main(config_f, start, end, chr):
    opts = get_opts(config_f)

    wd = opts["gen"]["wd"]
    run = Run(wd)

    exp1 = Experiment(opts["gen"]["f1"], opts["gen"]["type1"])
    exp2 = Experiment(opts["gen"]["f2"], opts["gen"]["type2"])

    plot = Plot(exp1, exp2, chr, start, end, run)
    print(plot.go())

    return 1

if __name__ == "__main__":
    config_f, start, end, chr = parse_args()
    sys.exit(main(config_f, start, end, chr))

