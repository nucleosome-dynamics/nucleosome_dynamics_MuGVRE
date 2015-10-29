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
from helpers import parse_exp_name, mkdir_p

###############################################################################

load_bams, nucler, nucdyn = ("{0}/{1}".format(RCODE, i)
                             for i in (LOAD_BAMS, NUCLER, NUCDYN))

###############################################################################

class Experiment:
    """
    Class to represent an MNase-seq experiment for nucleosome positioning
    """
    def __init__(self, fname, type):
        # experiments can be single-end or paired-end
        self.type = type
        # path to the experiment bam file
        self.bamfile = fname
        self.expname = parse_exp_name(fname)  # name of the experiment
        # where to store the intermediate RData with the preprocessed data
        self.rdatafile = self.expname + ".RData"
        # output file for nucleR
        self.nucler_out = self.expname + ".gff"

###############################################################################

def read_info_file(f):
    try:
        with open(f) as fh:
            done = set(line.strip() for line in fh)
    except FileNotFoundError:
        done = set()
    return done


def write_info_file(f, xs):
    with open(f, 'w') as fh:
        for x in xs:
            fh.write(x + '\n')

class Run:
    """
    Class to wrap the current estate of the pipeline and its calculations
    """
    def __init__(self, wd, done=None):
        self.wd = wd
        # create the working directory where files will be stored
        mkdir_p(self.wd)
        self.info_file = wd + "/info.txt"

        if done:
            self.done_calcs = done
            write_info_file(self.info_file, self.done_calcs)
        else:
            self.done_calcs = read_info_file(self.info_file)

    def write(self):
        """
        Write the current state of `self.done_calcs` into `info.txt`
        """
        write_info_file(self.info_file, self.done_calcs)

    def read(self):
        """
        Read the current state of `info.txt` into `self.done_calcs`
        """
        self.done_calcs = read_info_file(self.info_file)

    def add(self, x):
        """
        Add a given calculation to the set of finished calculations
        """
        self.done_calcs.add(x)
        self.write()

    def is_done(self, x):
        """
        Check if a calculation has finished
        """
        return x in self.done_calcs

###############################################################################

class Calculation:
    """
    Class to wrap the different possible kinds of calculations and
    ensure that things are only run once their calculation dependencies
    are met
    """
    def __init__(self, action, deps, run, calc_type, id=1):
        self.action = action
        self.deps = deps
        self.name = "{}{}".format(calc_type, id)
        self.run = run

    def go(self):
        if self.check_deps():
            self.action()
            self.run.add(self.name)
        else:
            print(("cannot run '{}' because at least one requiered " +
                   "step have not been performed").format(self.name))
            print("required steps:")
            for d in self.deps:
                print("\t" + d)
            sys.exit(1)

    def check_deps(self):
        return all(self.run.is_done(x) for x in self.deps)


class Load(Calculation):
    """
    Load the experiment into an RData file.
    """
    def __init__(self, exp, id, run):

        def load(exp, wd):
            subprocess.call([RUN_R, load_bams,
                             "--type", exp.type,
                             "--input", exp.bamfile,
                             "--output", "{}/{}".format(wd, exp.rdatafile)])

        self.deps = []
        self.name = "preproc{}".format(id)
        self.run = run
        self.action = lambda: load(exp, self.run.wd)


class NucleR(Calculation):
    """
    Run nucleR on an experiment
    """
    def __init__(self, exp, optargs, cores, id, run):

        def nucleR(exp, cores, wd, *args):
            subprocess.call([RUN_R, nucler,
                             "--input", "{}/{}".format(wd, exp.rdatafile),
                             "--cores", cores,
                             "--output", "{}/{}".format(wd, exp.nucler_out),
                             "--type", exp.type] + list(args))

        self.deps = ["preproc{}".format(id)]
        self.name = "nucleR{}".format(id)
        self.run = run
        self.action = lambda: nucleR(exp, cores, self.run.wd, *optargs)


class NucDyn(Calculation):
    """
    Run NucDyn on two Experiment objects
    """
    def __init__(self, exp1, exp2, optargs, cores, run):

        def nucleosome_dynamics(exp1, exp2, wd, cores, *args):
            fout = "{}/{}_{}.gff".format(wd, exp1.expname, exp2.expname)
            cmd = [RUN_R, nucdyn,
                   "--input1", "{}/{}".format(wd, exp1.rdatafile),
                   "--input2", "{}/{}".format(wd, exp2.rdatafile),
                   "--output", fout,
                   "--cores", cores] + list(args)
            subprocess.call(cmd)

        self.deps = ["preproc1", "preproc2"]
        self.name = "nucdyn"
        self.run = run
        self.action = lambda: nucleosome_dynamics(exp1, exp2,
                                                  self.run.wd, cores,
                                                  *optargs)

###############################################################################
