#!/usr/bin/env python3

import argparse
import itertools
import re
import sys


###############################################################################


class Nucleosome:
    """Class to represent a nucleosome call row"""
    def __init__(self, raw_str):
        self.raw_str = raw_str
        line_vals = raw_str.strip().split("\t")
        line_chr, _, _, line_start, line_end, line_score, *_ = line_vals
        self.chr = line_chr
        self.start = int(line_start)
        self.end = int(line_end)
        self.score = float(line_score)

    def in_range(self, start, end):
        none_check = start is None or end is None
        range_check = none_check or self.start >= start and self.end <= end
        return range_check

    def in_chr(self, chr):
        return chr is None or self.chr == chr

    def in_region(self, chr, start, end):
        none_chr = chr is None
        region_check = self.in_chr(chr) and self.in_range(start, end)
        return none_chr or region_check

    def __str__(self):
        return self.raw_str


def get_nucs(fh, chr=None, start=None, end=None):
    for line in fh:
        n = Nucleosome(line)
        if n.in_region(chr, start, end):
            yield n


def parse_range(x):
    pattern1 = r"(?P<chr>.+):(?P<start>\d+)\.\.(?P<end>\d+)$"
    pattern2 = r"^(?P<chr>[^:]+):?(.+)?$"
    if x:
        m = re.compile(pattern1).match(x)
        if m is not None:
            chr = m.group("chr")
            start = int(m.group("start"))
            end = int(m.group("end"))
        else:
            chr = re.compile(pattern2).match(x).group("chr")
            start = None
            end = None
    else:
        chr = None
        start = None
        end = None
    return chr, start, end


def format_range(chr=None, start=None, end=None):
    if chr:
        if start and end:
            return "{0}:{1}..{2}".format(chr, start, end)
        else:
            return chr
    else:
        return "ALL"


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("calls", help="gff file with the nucleR calls")
    parser.add_argument("output", help="output gff file")
    parser.add_argument("--range",
                        help="genomic range to process (ex: chrI:2000..2600)",
                        default="")
    args = parser.parse_args()
    return args.calls, args.output, args.range


###############################################################################


def mapcat(f, x):
    return list(itertools.chain.from_iterable(map(f, x)))


def rm_least(xs):
    i, _ = min(enumerate(x.score for x in xs), key=lambda a: a[1])
    xs.pop(i)
    return xs


def merged_rans(nucs):
    starts = [x.start for x in nucs]
    ends = [x.end for x in nucs]

    ii = [1] * len(starts)
    jj = [-1] * len(ends)

    stends = list(zip(starts+ends, ii+jj))
    stends.sort(key=lambda x: x[0])

    xs = list(itertools.accumulate(x for _, x in stends))

    si = (i for i, _ in enumerate(xs) if xs[i] != 0 and xs[i-1] == 0)
    ei = (i for i, x in enumerate(xs) if x == 0)

    cstarts = (stends[i][0] for i in si)
    cends = (stends[i][0] for i in ei)

    return list(zip(cstarts, cends))


def nuc_grouper(nucs):
    mrans = merged_rans(nucs)
    res = [[x for x in nucs if x.in_range(start, end)]
            for start, end in mrans]
    return res


def chr_unoverlapper(nucs):
    def do(x):
        if len(x) == 1:
            return x
        else:
            return mapcat(do, nuc_grouper(rm_least(x)))
    return mapcat(do, nuc_grouper(nucs))


def unoverlapper(nucs):
    def do(chr):
        return chr_unoverlapper([x for x in nucs if x.in_chr(chr)])
    chrs = set(x.chr for x in nucs)
    return mapcat(do, chrs)


###############################################################################


def main():
    calls_file, out_file, raw_range = get_args()

    ran_chr, ran_start, ran_end = parse_range(raw_range)

    print("input file:", calls_file)
    print("output file:", out_file)
    print("genomic range:", format_range(ran_chr, ran_start, ran_end))

    with open(calls_file) as fh:
        nucs = list(get_nucs(fh, ran_chr, ran_start, ran_end))

    res = unoverlapper(nucs)

    with open("/orozco/homes/pluto/rilla/test.gff", 'w') as fh:
        for nuc in res:
            fh.write(str(nuc))

    return 0


###############################################################################


if __name__ == "__main__":
    sys.exit(main())

###############################################################################
