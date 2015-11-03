#!/usr/bin/python3.4

"""
Definitions of some classes to represent different calculations
"""

###############################################################################

import subprocess
import sys

from defaults import RUN_R, RCODE, LOAD_BAMS, NUCLER, NUCDYN, PLOT

###############################################################################

def load(exp, wd):
    output = "{}/{}".format(wd, exp.rdatafile)
    subprocess.call([RUN_R, "{}/{}".format(RCODE, LOAD_BAMS),
                     "--type", exp.type,
                     "--input", exp.bamfile,
                     "--output", output])
    return output

def nucleR(exp, cores, wd, *args):
    output = "{}/{}".format(wd, exp.nucler_out)
    subprocess.call([RUN_R, "{}/{}".format(RCODE, NUCLER),
                     "--input", "{}/{}".format(wd, exp.rdatafile),
                     "--cores", cores,
                     "--output", output,
                     "--type", exp.type] + list(args))
    return output

def nucleosome_dynamics(exp1, exp2, wd, cores, *args):
    out_name = "{}/{}_{}".format(wd, exp1.expname, exp2.expname)
    fout_gff = out_name + ".gff"
    fout_rdata = out_name + ".RData"
    cmd = [RUN_R, "{}/{}".format(RCODE, NUCDYN),
           "--input1", "{}/{}".format(wd, exp1.rdatafile),
           "--input2", "{}/{}".format(wd, exp2.rdatafile),
           "--outputGff", fout_gff,
           "--outputRData", fout_rdata,
           "--cores", cores] + list(args)
    subprocess.call(cmd)
    return fout_gff, fout_rdata

def plot_dyn(exp1, exp2, wd, chr, start, end):
    base = "{}/{}_{}".format(wd, exp1.expname, exp2.expname)
    in_name = "{}.RData".format(base)
    out_name = "{}_{}_{}-{}.png".format(base, chr, start, end)
    cmd = [RUN_R, "{}/{}".format(RCODE, PLOT),
           "--input", in_name,
           "--output", out_name,
           "--start", str(start),
           "--end", str(end),
           "--chr", chr]
    subprocess.call(cmd)
    return out_name


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
            out = self.action()
            self.run.add(self.name)
            return out
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
        self.deps = []
        self.name = "preproc{}".format(id)
        self.run = run
        self.action = lambda: load(exp, self.run.wd)


class NucleR(Calculation):
    """
    Run nucleR on an experiment
    """
    def __init__(self, exp, optargs, cores, id, run):
        self.deps = ["preproc{}".format(id)]
        self.name = "nucleR{}".format(id)
        self.run = run
        self.action = lambda: nucleR(exp, cores, self.run.wd, *optargs)


class NucDyn(Calculation):
    """
    Run NucDyn on two Experiment objects
    """
    def __init__(self, exp1, exp2, optargs, cores, run):
        self.deps = ["preproc1", "preproc2"]
        self.name = "nucdyn"
        self.run = run
        self.action = lambda: nucleosome_dynamics(exp1, exp2,
                                                  self.run.wd, cores,
                                                  *optargs)


class Plot(Calculation):
    """
    Plot NucDyn for a region of the genome
    """
    def __init__(self, exp1, exp2, chr, start, end, run):
        self.deps = ["nucdyn"]
        self.run = run
        self.name = "plot_{}_{}-{}".format(chr, start, end)
        self.action = lambda: plot_dyn(exp1, exp2,
                                       self.run.wd,
                                       chr, start, end)

###############################################################################
