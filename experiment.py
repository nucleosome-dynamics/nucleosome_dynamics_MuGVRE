#!/usr/bin/python3.4

"""
Class that wraps the information of an experiment
"""

###############################################################################

class Experiment:
    """
    Class to represent an MNase-seq experiment for nucleosome positioning
    """
    def __init__(self, fname, type):

        def parse_exp_name(x):
            return x.split('/')[-1].split('.')[0]

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
