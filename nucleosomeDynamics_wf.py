#!/usr/bin/python3

import json
import os
import sys
import tarfile

from argparse   import ArgumentParser
from copy       import deepcopy
from functools  import partial
from itertools  import chain
from subprocess import call

###############################################################################

RPATH = "/usr/bin/Rscript"
BIN_BASE = "/orozco/services/Rdata/Web/apps/nucleServ_MuG"
PUBLIC_DIR = "/orozco/services/MuG/MuG_public/"

###############################################################################

def const(x, *_):
    return x

def dummy(*_):
    return []

###############################################################################

class StatsProc:
    @staticmethod
    def merger(*fhs):
        splitted = [[x.strip().split(',') for x in xs] for xs in fhs]
        ids = [x[0] for x in splitted[0]]
        content = [[x[1: ] for x in xs] for xs in splitted]
        def f(id, *xs):
            return ','.join(chain([id], chain.from_iterable(xs)))
        return '\n'.join(map(f, ids, *content))

    @staticmethod
    def gw_cleaner(fname, discard=('0', 'NA')):
        tmp_file = fname + "_tmp"
        with open(fname) as in_fh, open(tmp_file, 'w') as out_fh:
            for line in in_fh:
                _, *xs = line.strip().split(',')
                if not all(x in discard for x in xs):
                    out_fh.write(line)
        os.rename(tmp_file, fname)
        return fname

    @staticmethod
    def clean_genome_wide(fs):
        for f in fs:
            if f.endswith("_genes_stats.csv"):
                StatsProc.gw_cleaner(f)
        return fs

    @staticmethod
    def merge_tabs(x, gene_stats, col_order, w_dir):
        out_f = "{0}/{1}_genes_stats.csv".format(w_dir, x)
        potential_files = ["{0}/{1}_{2}_genes_stats.csv".format(w_dir, col, x)
                           for col in col_order]
        files_to_merge = [f for f in potential_files if f in gene_stats]

        if files_to_merge:
            in_fhs = map(open, files_to_merge)
            with open(out_f, 'w') as fh:
                fh.writelines(StatsProc.merger(*in_fhs))
            for fh in in_fhs:
                fh.close()
            return out_f, files_to_merge
        else:
            return None, []

    @staticmethod
    def merge_stats(stats_files, col_order):
        stats_files = deepcopy(stats_files)
        gene_stats = [x for x in stats_files if x.endswith("_genes_stats.csv")]
        tab_files = set("_".join(os.path.basename(x).split("_")[1:-2])
                        for x in gene_stats)
        w_dir = os.path.dirname(list(stats_files)[0])

        for x in tab_files:
            merged, to_merge = StatsProc.merge_tabs(x, gene_stats, col_order, w_dir)
            if merged is not None:
                stats_files.add(merged)
            for f in to_merge:
                stats_files.remove(f)
                os.remove(f)

        return stats_files

    @staticmethod
    def compress_stats(stat_files, in_files, out_dir):
        """
        Compress all the statistics files into a tgz file and return its
        metadata
        """
        output = os.path.join(out_dir, "statistics.tgz")

        with tarfile.open(output, 'w:gz') as fh:
            for f in stat_files:
                try:
                    fh.add(f, arcname=os.path.basename(f))
                    os.remove(f)
                except FileNotFoundError:
                    pass

        sources = list(set(x["value"] for x in in_files))
        return [{"name":      "statistics",
                 "file_path": output,
                 "sources":   sources}]

    @staticmethod
    def proc(stat_files, in_files, out_dir, col_order):
        merged_stats = StatsProc.merge_stats(stat_files, col_order)
        StatsProc.clean_genome_wide(merged_stats)
        stats_meta = StatsProc.compress_stats(merged_stats, in_files, out_dir)
        return stats_meta

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
        dir, fname = os.path.split(x)
        base, _ = os.path.splitext(fname)
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
        sources = [x["_id"] for x in xs]
        new_meta = deepcopy(meta)
        for m in new_meta:
            m["taxon_id"] = taxon_id
            m["sources"] = sources
            try:
                m["meta_data"]["assembly"] = assembly
            except KeyError:
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

class Calculation():
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
            os.makedirs(dir)
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
        call(cmd)
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
        return list(chain.from_iterable(res))


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
    """
    Convert the input BAM files into a temporary RData
    """
    exec_name = "readBAM"
    names = "MNaseSeq", "condition1", "condition2"

    def fun(self, f, public_dir, out_dir):
        input = f["file_path"]
        type = f["meta_data"]["paired"]

        in_dir, base = PathHelpers.base_name(input)
        output  = PathHelpers.build_path(base, in_dir,  "RData")

        args = {"input": input, "output": output, "type": type}
        meta = []
        return args, meta


@as_function
class nucleR(IterOnInfs, Bin):
    """
    Run nucleR for nucleosome positioning
    """
    exec_name = "nucleR"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        in_dir, base = PathHelpers.base_name(f["file_path"])
        input = PathHelpers.build_path(base, in_dir,  "RData")
        output = PathHelpers.build_path(base, out_dir, "gff", "NR")
        type = f["meta_data"]["paired"]

        args = {"input": input, "output": output, "type": type}
        meta = [{"name": "NR_gff", "file_path": output}]
        return args, meta


@as_function
class nfr(IterOnInfs, Bin):
    """
    Look for nucleosome-free regions
    """
    exec_name = "NFR"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = PathHelpers.base_name(f["file_path"])
        input = PathHelpers.build_path(base, out_dir, "gff", "NR")
        output = PathHelpers.build_path(base, out_dir, "gff", "NFR")

        args = {"input": input, "output": output}
        meta = [{"name": "NFR_gff", "file_path": output}]
        return args, meta


@as_function
class tss(IterOnInfs, Bin):
    """
    Classify the transcription start sites according to the nucleosomes that
    surround it
    """
    exec_name = "txstart"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        in_dir, base = PathHelpers.base_name(f["file_path"])

        calls = PathHelpers.build_path(base, out_dir, "gff", "NR")
        output = PathHelpers.build_path(base, out_dir, "gff", "TSS")

        assembly = f["meta_data"]["assembly"]
        genome = PathHelpers.get_genes_f(assembly, public_dir)

        args = {"calls": calls, "genome": genome, "output": output}
        meta = [{"name": "TSS_gff", "file_path": output}]
        return args, meta


@as_function
class period(IterOnInfs, Bin):
    """
    Calculate the periodicity and phase of the coverages on gene bodies
    """
    exec_name = "periodicity"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        in_dir, base = PathHelpers.base_name(f["file_path"])

        calls = PathHelpers.build_path(base, out_dir, "gff", "NR")
        reads = PathHelpers.build_path(base, in_dir,  "RData")
        gffOutput = PathHelpers.build_path(base, out_dir, "gff", "P")
        bwOutput = PathHelpers.build_path(base, out_dir, "bw", "P")

        assembly = f["meta_data"]["assembly"]
        type = f["meta_data"]["paired"]
        genes = PathHelpers.get_genes_f(assembly, public_dir)
        chrom_sizes = PathHelpers.get_chrom_sizes_f(assembly, public_dir)

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
    """
    Gaussian fittness and stiffness constant estimation
    """
    exec_name = "gausfitting"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        in_dir, base = PathHelpers.base_name(f["file_path"])

        reads = PathHelpers.build_path(base, in_dir, "RData")
        calls = PathHelpers.build_path(base, out_dir, "gff", "NR")
        output = PathHelpers.build_path(base, out_dir, "gff", "STF")

        args = {"calls": calls, "reads": reads, "output": output}
        meta = [{"name": "STF_gff", "file_path": output}]
        return args, meta


@as_function
class nuc_dyn(SelectTwo, Bin):
    """
    Compare two states with nucleosomeDynamics
    """
    exec_name = "nucDyn"
    name1 = "condition1"
    name2 = "condition2"

    def fun(self, f1, f2, public_dir, out_dir):
        splt = (PathHelpers.base_name(f["file_path"]) for f in (f1, f2))
        (in_dir1, base1), (in_dir2, base2) = splt
        nd_base = "{0}_{1}".format(base1, base2)

        input1, input2 = map(PathHelpers.build_path,
                             [base1, base2],
                             [in_dir1, in_dir2],
                             ["RData", "RData"])
        outputGff = PathHelpers.build_path(nd_base, out_dir, "gff", "ND")
        plotRData = PathHelpers.build_path(nd_base, out_dir, "RData", "ND", "plot")
        outputBigWig = PathHelpers.build_path(nd_base, out_dir, "bw", "ND")

        assembly = f1["meta_data"]["assembly"]
        genome = PathHelpers.get_chrom_sizes_f(assembly, public_dir)

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
    """
    nucleR's statistics
    """
    exec_name = "nucleR"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = PathHelpers.base_name(f["file_path"])

        input = PathHelpers.build_path(base, out_dir, "gff", "NR")
        out_genes = PathHelpers.build_path(base, out_dir, "csv", "NR", "genes_stats")
        out_gw = PathHelpers.build_path(base, out_dir, "csv", "NR", "stats")

        assembly = f["meta_data"]["assembly"]
        genome = PathHelpers.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "out_genes": out_genes,
                "out_gw":    out_gw,
                "genome":    genome}
        meta = [out_genes, out_gw]
        return args, meta


@as_function
class nfr_stats(IterOnInfs, Stats):
    """
    Nucleosome-free regions statistics
    """
    exec_name = "NFR"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = PathHelpers.base_name(f["file_path"])

        input = PathHelpers.build_path(base, out_dir, "gff", "NFR")
        out_gw = PathHelpers.build_path(base, out_dir, "csv", "NFR", "stats")

        assembly = f["meta_data"]["assembly"]
        genome = PathHelpers.get_genes_f(assembly, public_dir)

        args = {"input": input, "out_gw": out_gw, "genome": genome}
        meta = [out_gw]
        return args, meta


@as_function
class tss_stats(IterOnInfs, Stats):
    """
    Transcription start site classification statistics
    """
    exec_name = "txstart"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = PathHelpers.base_name(f["file_path"])

        input = PathHelpers.build_path(base, out_dir, "gff", "TSS")
        out_genes = PathHelpers.build_path(base, out_dir, "csv", "TSS", "genes_stats")
        out_gw = PathHelpers.build_path(base, out_dir, "png", "TSS", "stats1")
        out_gw2 = PathHelpers.build_path(base, out_dir, "png", "TSS", "stats2")

        assembly = f["meta_data"]["assembly"]
        genome = PathHelpers.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "genome":    genome,
                "out_genes": out_genes,
                "out_gw":    out_gw,
                "out_gw2":   out_gw2}
        meta = [out_genes, out_gw, out_gw2]
        return args, meta


@as_function
class period_stats(IterOnInfs, Stats):
    """
    Periodicity statistics
    """
    exec_name = "periodicity"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = PathHelpers.base_name(f["file_path"])

        input = PathHelpers.build_path(base, out_dir, "gff", "P")
        out_genes = PathHelpers.build_path(base, out_dir, "csv", "P", "genes_stats")
        out_gw = PathHelpers.build_path(base, out_dir, "csv", "P", "stats")

        assembly = f["meta_data"]["assembly"]
        genome = PathHelpers.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "genome":    genome,
                "out_genes": out_genes,
                "out_gw":    out_gw}
        meta = [out_genes, out_gw]
        return args, meta


@as_function
class gauss_stats(IterOnInfs, Stats):
    """
    Gaussian fitting and stiffness constant estimation statistics
    """
    exec_name = "gausfitting"
    names = "MNaseSeq",

    def fun(self, f, public_dir, out_dir):
        _, base = PathHelpers.base_name(f["file_path"])

        input = PathHelpers.build_path(base, out_dir, "gff", "STF")
        out_genes = PathHelpers.build_path(base, out_dir, "csv", "STF", "genes_stats")
        out_gw = PathHelpers.build_path(base, out_dir, "csv", "STF", "stats1")
        out_gw2 = PathHelpers.build_path(base, out_dir, "png", "STF", "stats2")

        assembly = f["meta_data"]["assembly"]
        genome = PathHelpers.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "genome":    genome,
                "out_genes": out_genes,
                "out_gw":    out_gw,
                "out_gw2":   out_gw2}
        meta = [out_genes, out_gw, out_gw2]
        return args, meta


@as_function
class nuc_dyn_stats(SelectTwo, Stats):
    """
    NucleosomeDynamics statistics
    """
    exec_name = "nucDyn"
    name1 = "condition1"
    name2 = "condition2"

    def fun(self, f1, f2, public_dir, out_dir):
        splt = (PathHelpers.base_name(f["file_path"]) for f in (f1, f2))
        (in_dir1, base1), (in_dir2, base2) = splt
        nd_base = "{0}_{1}".format(base1, base2)

        input = PathHelpers.build_path(nd_base, out_dir, "gff", "ND")
        out_genes = PathHelpers.build_path(nd_base, out_dir, "csv", "ND", "genes_stats")
        out_gw = PathHelpers.build_path(nd_base, out_dir, "png", "ND", "stats")

        assembly = f1["meta_data"]["assembly"]
        genome = PathHelpers.get_genes_f(assembly, public_dir)

        args = {"input":     input,
                "genome":    genome,
                "out_genes": out_genes,
                "out_gw":    out_gw}
        meta = [out_genes, out_gw]
        return args, meta

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
    parser.add_argument("--in_metadata",
                        required=True,
                        metavar="METADATA_JSON",
                        help="JSON file containing MuG metadata files")
    parser.add_argument("--out_metadata",
                        required=True,
                        metavar="RESULTS_JSON",
                        help="JSON file containing results metadata")
    parser.add_argument("--log_file",
                        required=False,
                        metavar="LOG_FILE",
                        help="Log file")


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


# def add_root(meta, root_dir):
    # """
    # Add the root directory to the path of each input file
    # """
    # d = dict(meta)
    # for v in d.values():
        # try:
            # v["file_path"] = os.path.join(root_dir, v["file_path"])
        # except KeyError:
            # pass
    # return d

###############################################################################

def cleanup(in_files, metadata):
    """
    Remove temporary files (i.e. RData files containing the reads)
    """
    for x in in_files:
        bam_file = metadata[x["value"]]["file_path"]
        base, _ = os.path.splitext(bam_file)
        rdata_file = "{0}.{1}".format(base, "RData")
        try:
            os.remove(rdata_file)
        except FileNotFoundError:
            pass

###############################################################################

class Calc:
    def __init__(self, name, bin=dummy, stats=dummy, deps=[], todo=False):
        self.name = name
        self.bin = bin
        self.stats = stats
        self.deps = deps
        self.todo = todo

    def mark_as_todo(self):
        self.todo = True
        for x in self.deps:
            x.mark_as_todo()

    def check_todo(self, xs):
        if self.name in xs:
            self.mark_as_todo()

    def run(self, in_files, metadata, arguments, public_dir, out_dir):
        calc_args = in_files, metadata, arguments, public_dir, out_dir
        if self.todo:
            return self.bin(*calc_args), self.stats(*calc_args)
        else:
            return [], []


class Run:
    def __init__(self, calcs, col_order):
        self.calcs = calcs
        self.col_order = col_order

    def run(self, in_files, metadata, arguments, public_dir, out_dir):
        asked = set(k for k, v in arguments.items() if v)
        calc_args = in_files, metadata, arguments, public_dir, out_dir
        for x in self.calcs:
            x.check_todo(asked)
        res = [x.run(*calc_args) for x in self.calcs]
        calcs_meta = list(chain.from_iterable(x for x, _ in res))
        stats_files = set(chain.from_iterable(x for _, x in res))
        stats_meta = StatsProc.proc(stats_files, in_files, out_dir, self.col_order)
        return {"output_files": list(chain(calcs_meta, stats_meta))}

###############################################################################

col_order = ("NR", "TSS", "P", "STF")

read_bam_calc = Calc("readBam",     read_bam)
nucleR_calc   = Calc("nucleR",      nucleR,  nucleR_stats,  [read_bam_calc])
nfr_calc      = Calc("NFR",         nfr,     nfr_stats,     [nucleR_calc])
tss_calc      = Calc("txstart",     tss,     tss_stats,     [nucleR_calc])
period_calc   = Calc("periodicity", period,  period_stats,  [nucleR_calc])
gauss_calc    = Calc("gausfitting", gauss,   gauss_stats,   [nucleR_calc])
nuc_dyn_calc  = Calc("nucDyn",      nuc_dyn, nuc_dyn_stats, [read_bam_calc])

my_calcs = (read_bam_calc,
            nucleR_calc,
            nfr_calc,
            tss_calc,
            period_calc,
            gauss_calc,
            nuc_dyn_calc)
my_run = Run(my_calcs, col_order)

###############################################################################

def main():
    # get arguments
    args = get_args()

    # load and parse inputs
    with open(args.config) as fh:
        config = json.load(fh)
    with open(args.in_metadata) as fh:
        in_meta = json.load(fh)

    out_metadata = args.out_metadata

    in_files = config["input_files"]
    arguments = get_args_dict(config["arguments"])

    metadata = preproc_meta(in_meta)
    out_dir = arguments["project"]

    out_meta = my_run.run(in_files, metadata, arguments, PUBLIC_DIR, out_dir)
    cleanup(in_files, metadata)

    json_out = json.dumps(out_meta, indent=4, separators=(',', ': '))
    with open(out_metadata, 'w') as fh:
        fh.write(json_out)

    return 0

###############################################################################

if __name__ == '__main__':
    sys.exit(main())

###############################################################################
