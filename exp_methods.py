#!/usr/bin/python3.4

"""
Definition of a class that represents an experiment and has a method to
perform a nucleosome calling on that experiment using nucleR.
Also contains a function to run NucDyn from two experiments
"""

###############################################################################

import os
import subprocess

from defaults import RUN_R, RCODE, LOAD_BAMS, NUCLER, NUCDYN, IN_DIR
from helpers import parse_exp_name

###############################################################################

load_bams, nucleR, nucdyn = ("{0}/{1}".format(RCODE, i)
                             for i in (LOAD_BAMS, NUCLER, NUCDYN))

###############################################################################

class Experiment:
    """
    Class to represent an MNase-seq experiment for nucleosome positioning
    """
    def __init__(self, fname, type, out_dir, loaded=False):
        # experiments can be single-end or paired-end
        self.type = type
        self.expname = parse_exp_name(fname)  # name of the experiment
        # path to the experiment bam file
        self.bamfile = "{0}/{1}.bam".format(IN_DIR, self.expname)
        # where to store the intermediate RData with the preprocessed data
        self.rdatafile = "{0}/{1}.RData".format(out_dir, self.expname)
        # output file for nucleR
        self.nucler_out = "{0}/{1}.gff".format(out_dir, self.expname)
        # flag to keep track of wheather the experiment has been loaded
        # and processed into an RData
        self.loaded = loaded

    def load(self):
        """
        Load the experiment into an RData file. Do nothing if it's already
        loaded.
        """
        if self.loaded:
            pass
        else:
            subprocess.call([RUN_R, load_bams,
                             "--type", self.type,
                             "--input", self.bamfile,
                             "--output", self.rdatafile])
            self.loaded = True

    def nucleR(self, cores, *args):
        """
        Run nucleR on that experiment
        """
        self.load()
        subprocess.call([RUN_R, nucleR,
                         "--input", self.rdatafile,
                         "--cores", cores,
                         "--output", self.nucler_out,
                         "--type", self.type] + list(args))

###############################################################################

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
