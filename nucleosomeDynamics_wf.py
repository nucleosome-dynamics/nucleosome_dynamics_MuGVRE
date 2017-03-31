#!/usr/bin/python3

import os
import sys
import argparse
import json
import re
import subprocess
import tarfile


#RPATH = "/opt/R-3.1.2/bin/Rscript"
RPATH = "/usr/bin/Rscript"
#BIN_BASE = "/home/rilla/nucleServ"
BIN_BASE = "/orozco/services/Rdata/Web/apps/nucleServ"


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


def get_bin_path(calc_name, calc_type="bin"):
    return "{0}/{1}/{2}.R".format(BIN_BASE, calc_type, calc_name)


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


def run_calc(calc, calc_type="bin", **kwargs):
    bin_path = get_bin_path(calc, calc_type=calc_type)
    arg_list = flatten([["--" + k, str(v)] for k, v in kwargs.items()])
    cmd = [RPATH, bin_path] + arg_list
    print("running", calc)
    subprocess.call(cmd)


def mkdir(dir):
    try:
        os.makedirs(dir)
    except FileExistsError:
        pass


def calculation(exec_name, calc_type="bin"):
    def decorator(fun):
        def g(*xs, args, public_dir, out_dir):
            new_args, out_meta = fun(*xs,
                                     public_dir=public_dir,
                                     out_dir=out_dir)
            if calc_type == "bin":
                calc_args = subset_args(exec_name, args)
                calc_args.update(new_args)
            elif calc_type == "statistics":
                calc_args = new_args
            mkdir(out_dir)
            run_calc(exec_name, calc_type=calc_type, **calc_args)
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


def build_path(base, root, extension, prefix=None, postfix=None):
    if prefix and postfix:
        return "{0}/{1}_{2}_{3}.{4}".format(root, prefix, base, postfix, extension)
    elif prefix:
        return "{0}/{1}_{2}.{3}".format(root, prefix, base, extension)
    elif postfix:
        return "{0}/{1}_{2}.{3}".format(root, base, postfix, extension)
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
@calculation("nucleR", "statistics")
def nucleR_stats(f, public_dir, out_dir):

    _, base = base_name(f["file_path"])
    input = build_path(base, out_dir, "gff", "NR")
    out_genes = build_path(base, out_dir, "csv", "NR", "genes_stats")
    out_gw = build_path(base, out_dir, "csv", "NR", "stats")
    assembly = f["meta_data"]["assembly"]
    genome = get_genes_f(assembly, public_dir)

    args = {"input":     input,
            "out_genes": out_genes,
            "out_gw":    out_gw,
            "genome":    genome}
    meta = [out_genes, out_gw]

    return args, meta


@iter_on_infs("MNaseSeq")
@calculation("NFR")
def nfr(f, public_dir, out_dir):

    _, base = base_name(f["file_path"])
    input  = build_path(base, out_dir, "gff", "NR")
    output = build_path(base, out_dir, "gff", "NFR")

    args = {"input": input, "output": output}
    meta = [{"name": "NFR_gff",
             "file_path": output,
             "source_id": [f["_id"]]}]

    return args, meta


@iter_on_infs("MNaseSeq")
@calculation("NFR", "statistics")
def nfr_stats(f, public_dir, out_dir):

    _, base = base_name(f["file_path"])
    input = build_path(base, out_dir, "gff", "NFR")
    out_gw = build_path(base, out_dir, "csv", "NFR", "stats")
    assembly = f["meta_data"]["assembly"]
    genome = get_genes_f(assembly, public_dir)

    args = {"input":     input,
            "out_gw":    out_gw,
            "genome":    genome}
    meta = [out_gw]

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
@calculation("txstart", "statistics")
def tss_stats(f, public_dir, out_dir):

    _, base = base_name(f["file_path"])
    input = build_path(base, out_dir, "gff", "TSS")
    out_genes = build_path(base, out_dir, "csv", "TSS", "genes_stats")
    out_gw = build_path(base, out_dir, "png", "TSS", "stats1")
    out_gw2 = build_path(base, out_dir, "png", "TSS", "stats2")
    assembly = f["meta_data"]["assembly"]
    genome = get_genes_f(assembly, public_dir)

    args = {"input":     input,
            "genome":    genome,
            "out_genes": out_genes,
            "out_gw":    out_gw,
            "out_gw2":   out_gw2}
    meta = [out_genes, out_gw, out_gw2]

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
@calculation("periodicity", "statistics")
def periodicity_stats(f, public_dir, out_dir):

    _, base = base_name(f["file_path"])
    input = build_path(base, out_dir, "gff", "P")
    out_genes = build_path(base, out_dir, "csv", "P", "genes_stats")
    out_gw = build_path(base, out_dir, "csv", "P", "stats")
    assembly = f["meta_data"]["assembly"]
    genome = get_genes_f(assembly, public_dir)

    args = {"input":     input,
            "genome":    genome,
            "out_genes": out_genes,
            "out_gw":    out_gw}
    meta = [out_genes, out_gw]

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


@iter_on_infs("MNaseSeq")
@calculation("gausfitting", "statistics")
def gauss_fit_stats(f, public_dir, out_dir):

    _, base = base_name(f["file_path"])
    input = build_path(base, out_dir, "gff", "STF")
    out_genes = build_path(base, out_dir, "csv", "STF", "genes_stats")
    out_gw = build_path(base, out_dir, "csv", "STF", "stats1")
    out_gw1 = build_path(base, out_dir, "png", "STF", "stats2")
    assembly = f["meta_data"]["assembly"]
    genome = get_genes_f(assembly, public_dir)

    args = {"input":     input,
            "genome":    genome,
            "out_genes": out_genes,
            "out_gw":    out_gw,
            "out_gw2":   out_gw2}
    meta = [out_genes, out_gw, out_gw2]

    return args, meta


@select_two("condition1", "condition2")
@calculation("nucDyn")
def nucDyn(f1, f2, public_dir, out_dir):

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


@iter_on_infs("MNaseSeq")
@calculation("nucDyn", "statistics")
def nucDyn_stats(f, public_dir, out_dir):

    _, base = base_name(f["file_path"])
    input = build_path(base, out_dir, "gff", "ND")
    out_genes = build_path(base, out_dir, "csv", "ND", "genes_stats")
    out_gw = build_path(base, out_dir, "png", "ND", "stats")
    assembly = f["meta_data"]["assembly"]
    genome = get_genes_f(assembly, public_dir)

    args = {"input":     input,
            "genome":    genome,
            "out_genes": out_genes,
            "out_gw":    out_gw}
    meta = [out_genes, out_gw]

    return args, meta


CALCS_ORDER = ["nucleR", "NFR", "txstart", "periodicity", "gausfitting", "nucDyn"]
CALCS = {"nucleR":      {"bin": nucleR,      "statistics": nucleR_stats},
         "NFR":         {"bin": nfr,         "statistics": nfr_stats},
         "txstart":     {"bin": tss,         "statistics": tss_stats},
         "periodicity": {"bin": periodicity, "statistics": periodicity_stats},
         "gausfitting": {"bin": gauss_fit,   "statistics": gauss_fit_stats},
         "nucDyn":      {"bin": nucDyn,      "statistics": nucDyn_stats}}


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


def make_stats(todo_calcs, in_files, metadata, arguments, public_dir, out_dir):

    calc_args = (in_files, metadata, arguments, public_dir, out_dir)
    stat_ls = (CALCS[i]["statistics"](*calc_args) for i in todo_calcs)
    stat_files = flatten(stat_ls)
    output = os.path.join(out_dir, "statistics.tar")

    if stat_files:
        with tarfile.open(output, 'w') as fh:
            for f in stat_files:
                try:
                    fh.add(f)
                    os.remove(f)
                except FileNotFoundError:
                    pass

        return [{"name": "statistics", "file_path": output}]
    else:
        return []


def cleanup(in_files, metadata):
    for x in in_files:
        bam_file = metadata[x["value"]]["file_path"]
        base, _ = os.path.splitext(bam_file)
        rdata_file = "{0}.{1}".format(base, "RData")
        try:
            os.remove(rdata_file)
        except FileNotFoundError:
            pass


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

    todo_calcs = find_calcs_todo(arguments)
    calc_args = (in_files, metadata, arguments, public_dir, out_dir)

    # run calculations
    calcs_meta = flatten(CALCS[i]["bin"](*calc_args) for i in todo_calcs)
    stats_meta = make_stats(todo_calcs, *calc_args)
    cleanup(in_files, metadata)

    # store output
    out_meta = calcs_meta + stats_meta
    json_out = json.dumps(out_meta, indent=4, separators=(',', ': '))
    with open(out_metadata, 'w') as fh:
        fh.write(json_out)

    return 0


if __name__ == '__main__':
    sys.exit(main())
