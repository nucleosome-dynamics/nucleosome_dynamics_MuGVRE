#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from functools import reduce
import operator as op
import subprocess
import sys

#good_chrs = {"chrI",    "chrII",  "chrIII",  "chrIV", "chrIX", "chrM", "chrV",
#             "chrVI",   "chrVII", "chrVIII", "chrX",  "chrXI", "chrXII",
#             "chrXIII", "chrXIV", "chrXV",   "chrXVI"}

SAMTOOLS = "/usr/bin/samtools"


def compose(*args):
    def comp2(f, g):
        return lambda x: f(g(x))
    return reduce(comp2, args)


def int2roman(x):
    assert type(x) == int
    assert 0 < x < 4000
    ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
    nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
    result = ""
    for i, n in zip(ints, nums):
        count = int(x / i)
        result += n * count
        x -= i * count
    return result


def add_prefix(chr):
    if chr.startswith("chr"):
        return chr
    else:
        return "chr" + chr


def fix_chr(chr):
    if chr in good_chrs:
        return chr


def make_roman(chr):
    try:
        return int2roman(int(chr))
    except ValueError:
        return chr


def check_mito(x):
    mito_options = {"mt", "mito", "m"}
    if x.lower() in mito_options:
        return "M"
    else:
        return x


def flatten(xss):
    return reduce(op.concat, xss)


def rm_prefix(chr):
    xs = ("chr", "chrom", "chromosome")
    prefixes = flatten(((x.capitalize(), x.upper(), x) for x in xs))
    fs = [(lambda x, p=p: x.replace(p, "")) for p in prefixes]
    return compose(*fs)(chr)


def fix_chr(chr):
    procs = (add_prefix,
             make_roman,
             check_mito,
             rm_prefix)
    return compose(*procs)(chr)


def proc_stream(stream):
    for line in stream:
        line = line.decode()
        if line.startswith("@SQ"):
            [chr] = [field.split(":")[1]
                     for field in line.split()
                     if field.startswith("SN")]
            yield line.replace(chr, fix_chr(chr))
        else:
            yield line


def get_args():
    descr = ("Try to fix names of the chromosomes in the header to make " +
             "them consistent")
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("-i",
                        "--input",
                        required=True,
                        help="bam input file")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        help="bam output file")
    args = parser.parse_args()
    return args.input, args.output


def rehead_bam(inf, outf):
    read_header = subprocess.Popen([SAMTOOLS, "view", "-H", inf],
                                   stdout=subprocess.PIPE)
    with open(outf, 'w') as fh:
        reheader = subprocess.Popen([SAMTOOLS, "reheader", "-", inf],
                                    stdin=subprocess.PIPE,
                                    stdout=fh)
        for line in proc_stream(read_header.stdout):
            reheader.stdin.write(line.encode())


def main():
    rehead_bam(*get_args())


if __name__ == "__main__":
    main()
