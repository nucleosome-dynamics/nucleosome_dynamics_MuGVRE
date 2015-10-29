#!/usr/bin/python3.4

"""
Definition of a class that represents an experiment and has a method to
perform a nucleosome calling on that experiment using nucleR.
Also contains a function to run NucDyn from two experiments
"""

###############################################################################

import os
import subprocess
import sys

from defaults import RUN_R, RCODE, LOAD_BAMS, NUCLER, NUCDYN
from helpers import parse_exp_name

###############################################################################

load_bams, nucleR, nucdyn = ("{0}/{1}".format(RCODE, i)
                             for i in (LOAD_BAMS, NUCLER, NUCDYN))

###############################################################################

class Experiment:
    """
    Class to represent an MNase-seq experiment for nucleosome positioning
    """
    def __init__(self, fname, type, out_dir):
        # experiments can be single-end or paired-end
        self.type = type
        # path to the experiment bam file
        self.bamfile = fname
        self.expname = parse_exp_name(fname)  # name of the experiment
        # where to store the intermediate RData with the preprocessed data
        self.rdatafile = "{0}/{1}.RData".format(out_dir, self.expname)
        # output file for nucleR
        self.nucler_out = "{0}/{1}.gff".format(out_dir, self.expname)


###############################################################################


def load(exp):
    """
    Load the experiment into an RData file. Do nothing if it's already
    loaded.
    """
    subprocess.call([RUN_R, load_bams,
                     "--type", exp.type,
                     "--input", exp.bamfile,
                     "--output", exp.rdatafile])


def nucleR(exp, cores, *args):
    """
    Run nucleR on that experiment
    """
    subprocess.call([RUN_R, nucleR,
                     "--input", exp.rdatafile,
                     "--cores", cores,
                     "--output", exp.nucler_out,
                     "--type", exp.type] + list(args))


def nucleosome_dynamics(exp1, exp2, out_dir, cores, *args):
    """
    Run NucDyn on two Experiment objects
    """
    fout = "{0}/{1}_{2}.gff".format(out_dir, exp1.expname, exp2.expname)
    exp1.load()
    exp2.load()
    cmd = [RUN_R, nucdyn,
           "--input1", exp1.rdatafile,
           "--input2", exp2.rdatafile,
           "--output", fout,
           "--cores", cores] + list(args)
    print(" ".join(cmd))
    subprocess.call(cmd)

###############################################################################

def read_info_file(f):
    try:
        with open(f) as fh:
            done = set(line.strip() for line in fh)
    except FileNotFoundError:
        done = set()
    return done

def write_info_file(f, x):
    done = read_info_file(f)
    if not x in done:
        with open(f, 'w') as fh:
            fh.write(x + "\n")

class Calculation:
    """
    Class to wrap the different possible kinds of calculations and
    ensure that things are only run once their calculation dependencies
    are met
    """
    def __init__(self, action, deps, log_f, calc_type, id=1):
        self.action = action
        self.deps = deps
        self.name = "{}{}".format(calc_type, id)
        self.log_f = log_f

    def run(self):
        if self.check_deps(self.log_f):
            self.action()
            write_info_file(self.log_f, self.name)
        else:
            print(("cannot run '{}' because at least one requiered " +
                   "step have not been performed").format(self.name))
            print("required steps:")
            for d in self.deps:
                print("\t" + d)
            sys.exit(1)

    def check_deps(self, f):
        done = read_info_file(self.log_f)
        return all((d in done for d in self.deps))


class Load(Calculation):
    def __init__(self, exp, id, log_f):
        self.action = lambda: load(exp)
        self.deps = []
        self.name = "load{}".format(id)
        self.log_f = log_f


class NucleR(Calculation):
    def __init__(self, exp, optargs, cores, id, log_f):
        self.action = lambda: nucleR(exp, cores, *optargs)
        self.deps = ["load{}".format(id)]
        self.name = "nucleR{}".format(id)
        self.log_f = log_f


class NucDyn(Calculation):
    def __init__(self, exp1, exp2, wd, optargs, cores, log_f):
        self.action = lambda: nucleosome_dynamics(exp1, exp2,
                                                  wd, cores,
                                                  *optargs)
        self.deps = ["load1", "load2"]
        self.name = "NucDyn"
        self.log_f = log_f
