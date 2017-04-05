#!/usr/bin/python3

from argparse import ArgumentParser
from collections import OrderedDict
from copy import deepcopy
from functools import partial
from itertools import chain
import json
from os import path, makedirs, remove
from subprocess import call
import sys
import tarfile

###############################################################################

#RPATH = "/opt/R-3.1.2/bin/Rscript"
RPATH = "/usr/bin/Rscript"
#BIN_BASE = "/home/rilla/nucleServ"
BIN_BASE = "/orozco/services/Rdata/Web/apps/nucleServ_MuG"

###############################################################################

class PathHelpers:
    """
    Class containing helpers for the calculations to format the input and
    output files
    """

    @staticmethod
    def get_chrom_sizes_f(assembly, public_dir):
        """
        Return the path of the chromosome sizes directory
        """
        return "{0}/refGenomes/{1}/{1}.fa.chrom.sizes".format(public_dir, assembly)

    @staticmethod
    def get_genes_f(assembly, public_dir):
        """
        Return the path of the gff containing the genes
        """
        return "{0}/refGenomes/{1}/genes.gff".format(public_dir, assembly)

    @staticmethod
    def base_name(x):
        """
        Given the path on an input file, return its directory and its base
        name without its extension
        """
        dir, fname = path.split(x)
        base, _ = path.splitext(fname)
        return dir, base

    @staticmethod
    def build_path(base, root, extension, prefix=None, postfix=None):
        """
        Given a base name, a directory, an extension and, optionally, a prefix
        and a postfix, construct the full name of an output file
        """
        if prefix and postfix:
            return "{0}/{1}_{2}_{3}.{4}".format(root, prefix, base, postfix, extension)
        elif prefix:
            return "{0}/{1}_{2}.{3}".format(root, prefix, base, extension)
        elif postfix:
            return "{0}/{1}_{2}.{3}".format(root, base, postfix, extension)
        else:
            return "{0}/{1}.{2}".format(root, base, extension)

###############################################################################

def const(x, *_):
    return x


class Bin:
    """
    Proper calculations (not statistics) should inherit from this class
    """
    type = "bin"

    def subset_args(self, args):
        """
        Given a full dictionary of arguments, get the ones specific to this
        calculation name, marked as `calc_name:arg_name`. Strip the prefix
        as well
        """
        prefix = self.exec_name + ":"
        res = {k.replace(prefix, ""): v
               for k, v in args.items()
               if k.startswith(prefix)}
        return res

    @staticmethod
    def add_meta(meta, *xs):
        """
        Add some metadata present in the input file to the metadata already
        returned by the calculation
        """
        x, *_ = xs
        assembly = x["meta_data"]["assembly"]
        taxon_id = x["taxon_id"]
        source_id = [x["_id"] for x in xs]
        new_meta = deepcopy(meta)
        for m in new_meta:
            m["taxon_id"] = taxon_id
            m["source_id"] = source_id
            try:
                m["meta_data"]["assembly"] = assembly
            except:
                m["meta_data"] = {"assembly": assembly}
        return new_meta

    def update_args(self, calc_args, args):
        """
        Add the arguments returned by the calculation preparing function
        (mostly input and output files) to the arguments already specified
        by the user
        """
        given_args = self.subset_args(args)
        new_args = deepcopy(calc_args)
        new_args.update(given_args)
        return new_args


class Stats:
    """
    Statistics run after calculations should inherit from this class
    """
    type = "statistics"
    add_meta = staticmethod(const)
    update_args = staticmethod(const)


class PreProc:
    """
    Used for data preprocessing (i.e. readBam)
    """
    type = "bin"
    add_meta = staticmethod(partial(const, []))
    update_args = staticmethod(const)

###############################################################################

class Calculation(PathHelpers):
    """
    All kinds of calculations should indirectly inherit from this (proper
    calculations and statistics). They should directly inherit from another
    class that deals with which inputs to compute.
    The calculations themselves should define their name, the input files they
    will use, and a function (`fun`) that returns the used command line
    arguments and metadata of its output files.
    """

    @staticmethod
    def mkdir(dir):
        """
        Create a directory if needed
        """
        try:
            makedirs(dir)
        except FileExistsError:
            pass

    def run_it(self, **kwargs):
        """
        Actually run the calculation with the supplied arguments
        """
        calc_type = self.type
        bin_path = self.get_bin_path(calc_type=calc_type)
        arg_pairs = (("--" + k, str(v)) for k, v in kwargs.items())
        arg_list = chain.from_iterable(arg_pairs)
        cmd = list(chain([RPATH, bin_path], arg_list))
        print_str = "running {0} ({1})".format(self.exec_name, calc_type)
        print(print_str)
        print(cmd)
        #call(cmd)
        print("==============================================================")

    def get_bin_path(self, calc_type="bin"):
        """
        Get the path of the script to run
        """
        return "{0}/{1}/{2}.R".format(BIN_BASE, calc_type, self.exec_name)

    def run_calc(self, *xs, args, public_dir, out_dir):
        """
        Run a calculation and return its output metadata
        """
        calc_args, out_meta = self.fun(*xs,
                                       public_dir=public_dir,
                                       out_dir=out_dir)
        out_meta = self.add_meta(out_meta, *xs)
        calc_args = self.update_args(calc_args, args)
        self.mkdir(out_dir)
        self.run_it(**calc_args)
        return out_meta

###############################################################################

class IterOnInfs(Calculation):
    """
    Given a list of inputs, the calculation will run on each one of them as
    long as they have the right name (a name present in the sequence
    self.names). Used for calculations that have only one input.
    """

    def run(self, in_files, metadata, args, public_dir, out_dir):
        ids = set([x["value"]
                  for x in in_files
                  if [x["name"] in self.names]])
        res = (self.run_calc(metadata[id],
                             args=args,
                             public_dir=public_dir,
                             out_dir=out_dir)
               for id in ids)
        return chain.from_iterable(res)


class SelectTwo(Calculation):
    """
    Given a list of inputs, the calculation will run on the first two that
    have the right names (corresponding to self.name1 and self.name2). Used for
    calculations that take two inputs (i.e. NucDyn).
    """

    @staticmethod
    def get_file(name, metadata, in_files):
        a, *_ = (metadata[i["value"]]
                 for i in in_files
                 if i["name"] == name)
        return a

    def run(self, in_files, metadata, args, public_dir, out_dir):
        f1, f2 = (self.get_file(n, metadata, in_files)
                  for n in (self.name1, self.name2))
        res = self.run_calc(f1,
                            f2,
                            args=args,
                            public_dir=public_dir,
                            out_dir=out_dir)
        return res

###############################################################################

def as_function(c):
    """
    Use this to decorate a class definition. With it, the class will define
    a function that returns the return value of its `run` method instead of
    returning an instance of an object
    """
    return c().run

###############################################################################

@as_function
class read_bam(IterOnInfs, PreProc):

    def __init__(self):
        self.exec_name = "readBAM"
        self.names = "MNaseSeq", "condition1", "condition2"

    def fun(self, f, public_dir, out_dir):
        input = f["file_path"]
        type = f["meta_data"]["paired"]

        in_dir, base = self.base_name(input)
        output  = self.build_path(base, in_dir,  "RData")

        args = {"input": input, "output": output, "type": type}
        meta = []
        return args, meta


@as_function
class nucleR(IterOnInfs, Bin):

    def __init__(self):
        self.exec_name = "nucleR"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        in_dir, base = self.base_name(f["file_path"])
        input = self.build_path(base, in_dir,  "RData")
        output = self.build_path(base, out_dir, "gff", "NR")
        type = f["meta_data"]["paired"]

        args = {"input": input, "output": output, "type": type}
        meta = [{"name": "NR_gff", "file_path": output}]
        return args, meta


@as_function
class nfr(IterOnInfs, Bin):

    def __init__(self):
        self.exec_name = "NFR"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = self.base_name(f["file_path"])
        input = self.build_path(base, out_dir, "gff", "NR")
        output = self.build_path(base, out_dir, "gff", "NFR")

        args = {"input": input, "output": output}
        meta = [{"name": "NFR_gff", "file_path": output}]
        return args, meta


@as_function
class tss(IterOnInfs, Bin):

    def __init__(self):
        self.exec_name = "txstart"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        in_dir, base = self.base_name(f["file_path"])

        calls = self.build_path(base, out_dir, "gff", "NR")
        output = self.build_path(base, out_dir, "gff", "TSS")

        assembly = f["meta_data"]["assembly"]
        genome = self.get_genes_f(assembly, public_dir)

        args = {"calls": calls, "genome": genome, "output": output}
        meta = [{"name": "TSS_gff", "file_path": output}]
        return args, meta


@as_function
class period(IterOnInfs, Bin):

    def __init__(self):
        self.exec_name = "periodicity"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        in_dir, base = self.base_name(f["file_path"])

        calls = self.build_path(base, out_dir, "gff", "NR")
        reads = self.build_path(base, in_dir,  "RData")
        gffOutput = self.build_path(base, out_dir, "gff", "P")
        bwOutput = self.build_path(base, out_dir, "bw", "P")

        assembly = f["meta_data"]["assembly"]
        type = f["meta_data"]["paired"]
        genes = self.get_genes_f(assembly, public_dir)
        chrom_sizes = self.get_chrom_sizes_f(assembly, public_dir)

        args = {"calls":       calls,
                "reads":       reads,
                "type":        type,
                "gffOutput":   gffOutput,
                "bwOutput":    bwOutput,
                "genes":       genes,
                "chrom_sizes": chrom_sizes}
        meta = [{"name": "P_gff", "file_path": gffOutput},
                {"name": "P_bw",  "file_path": bwOutput}]
        return args, meta


@as_function
class gauss(IterOnInfs, Bin):

    def __init__(self):
        self.exec_name = "gausfitting"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        in_dir, base = self.base_name(f["file_path"])

        reads = self.build_path(base, in_dir, "RData")
        calls = self.build_path(base, out_dir, "gff", "NR")
        output = self.build_path(base, out_dir, "gff", "STF")

        args = {"calls": calls, "reads": reads, "output": output}
        meta = [{"name": "STF_gff", "file_path": output}]
        return args, meta


@as_function
class nuc_dyn(SelectTwo, Bin):

    def __init__(self):
        self.exec_name = "nucDyn"
        self.name1 = "condition1"
        self.name2 = "condition2"

    def fun(self, f1, f2, public_dir, out_dir):
        splt = (self.base_name(f["file_path"]) for f in (f1, f2))
        (in_dir1, base1), (in_dir2, base2) = splt
        nd_base = "{0}_{1}".format(base1, base2)

        input1, input2 = map(self.build_path,
                             [base1, base2],
                             [in_dir1, in_dir2],
                             ["RData", "RData"])
        outputGff = self.build_path(nd_base, out_dir, "gff", "ND")
        plotRData = self.build_path(nd_base, out_dir, "RData", "ND")
        outputBigWig = self.build_path(nd_base, out_dir, "bw", "ND")

        assembly = f1["meta_data"]["assembly"]
        genome = self.get_chrom_sizes_f(assembly, public_dir)

        args = {"input1":       input1,
                "input2":       input2,
                "outputGff":    outputGff,
                "outputBigWig": outputBigWig,
                "plotRData":    plotRData,
                "genome":       genome}
        meta = [{"name": "ND_gff", "file_path": outputGff},
                {"name": "ND_bw",  "file_path": outputBigWig}]
        return args, meta

###############################################################################

@as_function
class nucleR_stats(IterOnInfs, Stats):

    def __init__(self):
        self.exec_name = "nucleR"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = self.base_name(f["file_path"])

        input = self.build_path(base, out_dir, "gff", "NR")
        out_genes = self.build_path(base, out_dir, "csv", "NR", "genes_stats")
        out_gw = self.build_path(base, out_dir, "csv", "NR", "stats")

        assembly = f["meta_data"]["assembly"]
        genome = self.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "out_genes": out_genes,
                "out_gw":    out_gw,
                "genome":    genome}
        meta = [out_genes, out_gw]
        return args, meta


@as_function
class nfr_stats(IterOnInfs, Stats):

    def __init__(self):
        self.exec_name = "NFR"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = self.base_name(f["file_path"])

        input = self.build_path(base, out_dir, "gff", "NFR")
        out_gw = self.build_path(base, out_dir, "csv", "NFR", "stats")

        assembly = f["meta_data"]["assembly"]
        genome = self.get_genes_f(assembly, public_dir)

        args = {"input": input, "out_gw": out_gw, "genome": genome}
        meta = [out_gw]
        return args, meta


@as_function
class tss_stats(IterOnInfs, Stats):

    def __init__(self):
        self.exec_name = "txstart"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = self.base_name(f["file_path"])

        input = self.build_path(base, out_dir, "gff", "TSS")
        out_genes = self.build_path(base, out_dir, "csv", "TSS", "genes_stats")
        out_gw = self.build_path(base, out_dir, "png", "TSS", "stats1")
        out_gw2 = self.build_path(base, out_dir, "png", "TSS", "stats2")

        assembly = f["meta_data"]["assembly"]
        genome = self.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "genome":    genome,
                "out_genes": out_genes,
                "out_gw":    out_gw,
                "out_gw2":   out_gw2}
        meta = [out_genes, out_gw, out_gw2]
        return args, meta


@as_function
class period_stats(IterOnInfs, Stats):

    def __init__(self):
        self.exec_name = "periodicity"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = self.base_name(f["file_path"])

        input = self.build_path(base, out_dir, "gff", "P")
        out_genes = self.build_path(base, out_dir, "csv", "P", "genes_stats")
        out_gw = self.build_path(base, out_dir, "csv", "P", "stats")

        assembly = f["meta_data"]["assembly"]
        genome = self.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "genome":    genome,
                "out_genes": out_genes,
                "out_gw":    out_gw}
        meta = [out_genes, out_gw]
        return args, meta


@as_function
class gauss_stats(IterOnInfs, Stats):

    def __init__(self):
        self.exec_name = "gausfitting"
        self.names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = self.base_name(f["file_path"])

        input = self.build_path(base, out_dir, "gff", "STF")
        out_genes = self.build_path(base, out_dir, "csv", "STF", "genes_stats")
        out_gw = self.build_path(base, out_dir, "csv", "STF", "stats1")
        out_gw2 = self.build_path(base, out_dir, "png", "STF", "stats2")

        assembly = f["meta_data"]["assembly"]
        genome = self.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "genome":    genome,
                "out_genes": out_genes,
                "out_gw":    out_gw,
                "out_gw2":   out_gw2}
        meta = [out_genes, out_gw, out_gw2]
        return args, meta


@as_function
class nuc_dyn_stats(SelectTwo, Stats):

    def __init__(self):
        self.exec_name = "nucDyn"
        self.name1 = "condition1"
        self.name2 = "condition2"

    def fun(self, f1, f2, public_dir, out_dir):
        splt = (self.base_name(f["file_path"]) for f in (f1, f2))
        (in_dir1, base1), (in_dir2, base2) = splt
        nd_base = "{0}_{1}".format(base1, base2)

        input = self.build_path(nd_base, out_dir, "gff", "ND")
        out_genes = self.build_path(nd_base, out_dir, "csv", "ND", "genes_stats")
        out_gw = self.build_path(nd_base, out_dir, "png", "ND", "stats")

        assembly = f1["meta_data"]["assembly"]
        genome = self.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "genome":    genome,
                "out_genes": out_genes,
                "out_gw":    out_gw}
        meta = [out_genes, out_gw]
        return args, meta

###############################################################################

def dummy(*_):
    return []

CALCS = OrderedDict(
    [("readBam",     {"bin": read_bam, "stats": dummy,         "deps": []}),
     ("nucleR",      {"bin": nucleR,   "stats": nucleR_stats,  "deps": ["readBam"]}),
     ("NFR",         {"bin": nfr,      "stats": nfr_stats,     "deps": ["nucleR"]}),
     ("txstart",     {"bin": tss,      "stats": tss_stats,     "deps": ["nucleR"]}),
     ("periodicity", {"bin": period,   "stats": period_stats,  "deps": ["nucleR"]}),
     ("gausfitting", {"bin": gauss,    "stats": gauss_stats,   "deps": ["nucleR"]}),
     ("nucDyn",      {"bin": nuc_dyn,  "stats": nuc_dyn_stats, "deps": ["readBam"]})
])

###############################################################################

def get_args():
    """
    Parse the command line arguments
    """
    parser = ArgumentParser(prog="nucleosomeDynamics_wf",
                            description="Nucleoseom Dynamics workflow")

    parser.add_argument("--config",
                        required=True,
                        metavar="CONFIG_JSON",
                        help="JSON file containing workflow parameters")
    parser.add_argument("--root_dir",
                        required=True,
                        metavar="ABS_PATH",
                        help="Absolute path of the user data directory.")
    parser.add_argument("--public_dir",
                        required=False,
                        metavar="PUBLIC_PATH",
                        help="Absolute path of the MuG public directory (with reference genome data, etc).")
    parser.add_argument("--metadata",
                        required=True,
                        metavar="METADATA_JSON",
                        help="JSON file containing MuG metadata files")
    parser.add_argument("--out_metadata",
                        required=False,
                        metavar="RESULTS_JSON",
                        help="JSON file containing results metadata")

    return parser.parse_args()

###############################################################################

def preproc_meta(metadata):
    """
    For easier access, convert the metadata list into a dictionary where
    the ids are the keys
    """
    res = {}
    for x in metadata:
        try:
            k = x["_id"]
        except KeyError:
            continue
        res[k] = x
    return res


def get_args_dict(xs):
    """
    For easier access, convert the format of the arguments from a list of
    dictionaries with the values `name` and `value` into a single dictionary
    where the values of `name` are the keys and the values of `value` are the
    values
    """
    return {x["name"]: x["value"] for x in xs}


def find_calcs_todo(args):
    """
    Return a two lists of the calculations to be run in order. The first list
    constains the functions for proper calculations and the second one is for
    the statistics
    """

    def try_get_val(d, i):
        try:
            return d[i]
        except KeyError:
            return False

    def set_asked_calcs(calc_ls, asked_ls):
        if asked_ls:
            for k, v in calc_ls.items():
                if k in asked_ls:
                    v["todo"] = True
                    set_asked_calcs(calc_ls, v["deps"])
            return calc_ls
        else:
            return calc_ls

    my_calcs = deepcopy(CALCS)
    for v in my_calcs.values():
        v["todo"] = False
    asked_calcs = [k for k in my_calcs.keys() if try_get_val(args, k)]
    set_asked_calcs(my_calcs, asked_calcs)
    bins = (x["bin"] for x in my_calcs.values() if x["todo"])
    stats = (x["stats"] for x in my_calcs.values() if x["todo"])
    return bins, stats


def add_root(meta, root_dir):
    """
    Add the root directory to the path of each input file
    """
    d = dict(meta)
    for v in d.values():
        try:
            v["file_path"] = path.join(root_dir, v["file_path"])
        except KeyError:
            pass
    return d

###############################################################################

def make_stats(stat_files, in_files, out_dir):
    """
    Compress all the statistics files into a tgz file and return its metadata
    """
    output = path.join(out_dir, "statistics.tgz")

    with tarfile.open(output, 'w:gz') as fh:
        for f in stat_files:
            try:
                fh.add(f, arcname=path.basename(f))
                remove(f)
            except FileNotFoundError:
                pass

    source_id = list(set(x["value"] for x in in_files))
    return [{"name":      "statistics",
             "file_path": output,
             "source_id": source_id}]


def cleanup(in_files, metadata):
    """
    Remove temporary files (i.e. RData files containing the reads)
    """
    for x in in_files:
        bam_file = metadata[x["value"]]["file_path"]
        base, _ = path.splitext(bam_file)
        rdata_file = "{0}.{1}".format(base, "RData")
        try:
            remove(rdata_file)
        except FileNotFoundError:
            pass

###############################################################################

def main():
    # get arguments
    args = get_args()

    # load and parse inputs
    with open(args.config) as fh:
        config = json.load(fh)
    with open(args.metadata) as fh:
        meta = json.load(fh)

    root_dir = args.root_dir
    public_dir = args.public_dir
    out_metadata = args.out_metadata

    in_files = config["input_files"]
    arguments = get_args_dict(config["arguments"])
    metadata = add_root(preproc_meta(meta), root_dir)
    out_dir = path.join(root_dir, arguments["project"])
    calc_args = (in_files, metadata, arguments, public_dir, out_dir)

    bin_funs, stats_funs = find_calcs_todo(arguments)
    # run calculations
    calcs_res = (f(*calc_args) for f in bin_funs)
    # we need to force evaluation here so that calcs run before stats
    calcs_meta = list(chain.from_iterable(calcs_res))
    stat_ls = (f(*calc_args) for f in stats_funs)
    stat_files = chain.from_iterable(stat_ls)
    stats_meta = make_stats(stat_files, in_files, out_dir)

    cleanup(in_files, metadata)

    # store output
    out_meta = {"output_files": list(chain(calcs_meta, stats_meta))}
    json_out = json.dumps(out_meta, indent=4, separators=(',', ': '))
    with open(out_metadata, 'w') as fh:
        fh.write(json_out)

    return 0

###############################################################################

if __name__ == '__main__':
    sys.exit(main())
