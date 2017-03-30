#!/usr/bin/python2.7

import os
import sys
import argparse
import json
import re
import subprocess


RPATH = "/opt/R-3.1.2/bin/Rscript"
BIN_BASE = "/home/rilla/nucleServ/bin/{0}.R"


def get_args():
    parser = argparse.ArgumentParser(prog="nucleosomeDynamics_wf",
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


def subset_args(x, args):
    prefix = x + ":"
    res = {k.replace(prefix, ""): v
           for k, v in args.items()
           if k.startswith(prefix)}
    return res


def preproc_meta(metadata):
    def f(xs):
        for x in xs:
            try:
                k = x["_id"]
                yield k, x
            except KeyError:
                continue
    return {k: v for k, v in f(metadata)}


def get_args_dict(xs):
    return {x["name"]: x["value"] for x in xs}


def flatten(xss):
    return [x for xs in xss for x in xs]


def get_genome_path(my_assembly, public_dir):
    base_dir = os.path.join(public_dir, "refGenomes")
    for organism in os.listdir(base_dir):
        organism_dir = os.path.join(base_dir, organism)
        for assembly in os.listdir(organism_dir):
            if my_assembly == assembly:
                return os.path.join(organism_dir, assembly)
    return None


def get_chrom_sizes_f(assembly, public_dir):
    genome_path = get_genome_path(assembly, public_dir)
    return "{0}/{1}.fa.chrom.sizes".format(genome_path, assembly)


def get_genes_f(assembly, public_dir):
    genome_path = get_genome_path(assembly, public_dir)
    return os.path.join(genome_path, "genes.gff")


def run_calc(calc, **kwargs):
    bin_path = BIN_BASE.format(calc)
    arg_list = flatten([["--" + k, str(v)] for k, v in kwargs.items()])
    cmd = [RPATH, bin_path] + arg_list
    #print(cmd)
    #subprocess.call(cmd)


def mkdir(dir):
    try:
        os.makedirs(dir)
    except FileExistsError:
        pass


def calculation(exec_name):
    def decorator(fun):
        def g(*xs, args, public_dir, out_dir):
            calc_args = subset_args(exec_name, args)
            new_args, out_meta = fun(*xs,
                                     public_dir=public_dir,
                                     out_dir=out_dir)
            calc_args.update(new_args)
            mkdir(out_dir)
            run_calc(exec_name, **calc_args)
            return out_meta
        return g
    return decorator


def iter_on_infs(*names):
    def decorator(fun):
        def g(in_files, metadata, args, public_dir, out_dir):
            res = [fun(metadata[x["value"]],
                       args=args,
                       public_dir=public_dir,
                       out_dir=out_dir)
                   for x in in_files
                   if x["name"] in names]
            return flatten(res)
        return g
    return decorator


def select_two(x, y):
    def decorator(fun):
        def g(in_files, metadata, args, public_dir, out_dir):
            def get_file(name):
                a, *_ = (metadata[i["value"]]
                         for i in in_files
                         if i["name"] == name)
                return a
            f1, f2 = (get_file(n) for n in (x, y))
            res = fun(f1,
                      f2,
                      args=args,
                      public_dir=public_dir,
                      out_dir=out_dir)
            return res
        return g
    return decorator


def base_name(x):
    dir, fname = os.path.split(x)
    base, _ = os.path.splitext(fname)
    return dir, base


def build_path(base, root, extension, prefix=None):
    if prefix:
        return "{0}/{1}_{2}.{3}".format(root, prefix, base, extension)
    else:
        return "{0}/{1}.{2}".format(root, base, extension)


def build_ext_paths(base, root, extensions, prefix):
    return (build_path(base, root, x, prefix) for x in extensions)


@iter_on_infs("MNaseSeq", "condition1", "condition2")
@calculation("readBAM")
def read_bam(f, public_dir, out_dir):

    input = f["file_path"]
    type = f["meta_data"]["paired"]
    in_dir, base = base_name(input)
    output  = build_path(base, in_dir,  "RData")

    args = {"input": input, "output": output, "type": type}
    meta = []

    return args, meta


@iter_on_infs("MNaseSeq")
@calculation("nucleR")
def nucleR(f, public_dir, out_dir):

    in_dir, base = base_name(f["file_path"])
    input  = build_path(base, in_dir,  "RData")
    output = build_path(base, out_dir, "gff", "NR")
    type = f["meta_data"]["paired"]
    assembly = f["meta_data"]["assembly"]

    args = {"input": input, "output": output, "type": type}
    meta = [{"name":      "NR_gff",
             "file_path": output,
             "source_id": [f["_id"]],
             "meta_data": {"assembly": assembly}}]

    return args, meta


@iter_on_infs("MNaseSeq")
@calculation("NFR")
def nfr(f, _, out_dir):

    _, base = base_name(f["file_path"])
    input  = build_path(base, out_dir, "gff", "NR")
    output = build_path(base, out_dir, "gff", "NFR")

    args = {"input": input, "output": output}
    meta = [{"name": "NFR_gff",
             "file_path": output,
             "source_id": [f["_id"]]}]

    return args, meta


@iter_on_infs("MNaseSeq")
@calculation("txstart")
def tss(f, public_dir, out_dir):

    in_dir, base = base_name(f["file_path"])
    calls    = build_path(base, out_dir, "gff", "NR")
    output   = build_path(base, out_dir, "gff", "TSS")
    assembly = f["meta_data"]["assembly"]
    genome = get_genes_f(assembly, public_dir)

    args = {"calls":    calls,
            "genome":   genome,
            "output":   output}
    meta = [{"name":     "TSS_gff",
             "file_path": output,
             "source_id": [f["_id"]]}]

    return args, meta


@iter_on_infs("MNaseSeq")
@calculation("periodicity")
def periodicity(f, public_dir, out_dir):

    in_dir, base = base_name(f["file_path"])
    calls = build_path(base, out_dir, "gff", "NR")
    reads = build_path(base, in_dir,  "RData")
    gffOutput, bwOutput = build_ext_paths(base, out_dir, ["gff", "bw"], "P")
    assembly = f["meta_data"]["assembly"]
    type = f["meta_data"]["paired"]
    genes = get_genes_f(assembly, public_dir)
    chrom_sizes = get_chrom_sizes_f(assembly, public_dir)

    args = {"calls":       calls,
            "reads":       reads,
            "type":        type,
            "gffOutput":   gffOutput,
            "bwOutput":    bwOutput,
            "genes":       genes,
            "chrom_sizes": chrom_sizes}

    meta = [{"name":      "P_gff",
             "file_path": gffOutput,
             "source_id": [f["_id"]]},
            {"name":      "P_bw",
             "file_path": bwOutput,
             "source_id": [f["_id"]]}]

    return args, meta


@iter_on_infs("MNaseSeq")
@calculation("gausfitting")
def gauss_fit(f, _, out_dir):

    in_dir, base = base_name(f["file_path"])
    reads = build_path(base, in_dir, "RData")
    calls = build_path(base, out_dir, "gff", "NR")
    output = build_path(base, out_dir, "gff", "STF")

    args = {"calls": calls, "reads": reads, "output": output}
    meta = [{"name":     "STF_gff",
             "file_path": output,
             "source_id": [f["_id"]]}]

    return args, meta


@select_two("condition1", "condition2")
@calculation("nucDyn")
def nucleosomeDynamics(f1, f2, public_dir, out_dir):

    splt = (base_name(f["file_path"]) for f in (f1, f2))
    (in_dir1, base1), (in_dir2, base2) = splt
    input1, input2 = map(build_path,
                         [base1, base2],
                         [in_dir1, in_dir2],
                         ["RData", "RData"])
    nd_base = "{0}_{1}".format(base1, base2)
    outs = build_ext_paths(nd_base, out_dir, ["gff", "RData", "bw"], "ND")
    outputGff, plotRData, outputBigWig = outs
    assembly = f1["meta_data"]["assembly"]
    genome = get_chrom_sizes_f(assembly, public_dir)

    args = {"input1":       input1,
            "input2":       input2,
            "outputGff":    outputGff,
            "outputBigWig": outputBigWig,
            "plotRData":    plotRData,
            "genome":       genome}
    meta = [{"name":      "ND_gff",
             "file_path": outputGff,
             "source_id": [f1["_id"], f2["_id"]]},
            {"name":      "ND_bw",
             "file_path": outputBigWig,
             "source_id": [f1["_id"], f2["_id"]]}]

    return args, meta


CALCS_ORDER = ["nucleR", "NFR", "txstart", "periodicity", "gausfitting", "nucDyn"]
CALCS = {"nucleR":      nucleR,
         "NFR":         nfr,
         "txstart":     tss,
         "periodicity": periodicity,
         "gausfitting": gauss_fit,
         "nucDyn":      nucleosomeDynamics}


def find_calcs_todo(args):

    def try_get_val(d, i):
        try:
            return d[i]
        except KeyError:
            return False

    return [x for x in CALCS_ORDER if try_get_val(args, x)]


def add_root(meta, root_dir):
    d = dict(meta)
    for v in d.values():
        try:
            v["file_path"] = os.path.join(root_dir, v["file_path"])
        except KeyError:
            pass
    return d


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
    out_dir = os.path.join(root_dir, arguments["project"])

    # load bam files
    read_bam(in_files, metadata, arguments, public_dir, out_dir)

    # run calculations
    todo_calcs = find_calcs_todo(arguments)
    out_meta = flatten(CALCS[i](in_files,
                                metadata,
                                arguments,
                                public_dir,
                                out_dir)
                       for i in todo_calcs)

    # cleanup
    for x in in_files:
        bam_file = metadata[x["value"]]["file_path"]
        base, _ = os.path.splitext(bam_file)
        rdata_file = "{0}.{1}".format(base, "RData")
        try:
            os.remove(rdata_file)
        except FileNotFoundError:
            pass

    # store output
    json_out = json.dumps(out_meta, indent=4, separators=(',', ': '))
    with open(out_metadata, 'w') as fh:
        fh.write(json_out)

    return 0


if __name__ == '__main__':
    sys.exit(main())
