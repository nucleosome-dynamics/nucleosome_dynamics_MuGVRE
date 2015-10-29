#!/usr/bin/python3.4

"""
Defition of a class that wraps the state of the whole pipeline
"""

###############################################################################

import os

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

def mkdir_p(d):
    try:
        os.makedirs(d)
    except FileExistsError:
        pass


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
