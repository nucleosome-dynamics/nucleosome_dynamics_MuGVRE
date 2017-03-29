#!/usr/bin/env python3

import argparse
import itertools
import re
import sys

###############################################################################

def parse_range(x, pattern=r"(?P<chr>.+):(?P<start>\d+)\.\.(?P<end>\d+)"):
    match = re.compile(pattern).match(x)
    chr = match.group("chr")
    start = int(match.group("start"))
    end = int(match.group("end"))
    return chr, start, end


def parse_line(line):
    line_chr, _, _, line_start, line_end, *_ = line.split("\t")
    return line_chr, int(line_start), int(line_end)


def filter_lines(x, chr, start, end):
    c, s, e, = x
    return chr == c and s >= start and e <= end


def get_nucs(fh, chr, start, end):
    for line in fh:
        parsed = parse_line(line)
        if filter_lines(parsed, chr, start, end):
            _, s, e = parsed
            yield s, e


def read_fasta(f):
    def iter(xs):
        for x in xs:
            lines = x.split()
            chr = lines[0]
            seq = "".join(lines[1:])
            yield chr, seq
    with open(f) as fh:
        txt = fh.read()
    res = {k: v for k, v in iter(txt.split(">")[1:])}
    return res


def get_seq(fh):
    return ''.join(line.strip() for line in fh)


def firsts(xs):
    return (x for x, *_ in xs)


def seconds(xs):
    return (x for _, x, *_ in xs)


def get_range(nucs):
    start = min(firsts(nucs))
    end = max(seconds(nucs))
    return start, end


def shift(x, by=0):
    start, end = x
    return start-by, end-by


def expand_range(x, by=0):
    start, end = x
    return start-by, end+by


def merge(ranges):
    ranges = list(ranges)
    starts = list(firsts(ranges))
    ends = list(seconds(ranges))

    ii = [1] * len(starts)
    jj = [-1] * len(ends)

    stends = list(zip(starts+ends, ii+jj))
    stends.sort(key=lambda x: x[0])

    xs = list(itertools.accumulate(seconds(stends)))

    si = (i for i, _ in enumerate(xs) if xs[i] != 0 and xs[i-1] == 0)
    ei = (i for i, x in enumerate(xs) if x == 0)

    cstarts = (stends[i][0] for i in si)
    cends = (stends[i][0] for i in ei)

    return list(zip(cstarts, cends))


def cut_seq(seq, ranges):
    ii = [0] + list(seconds(ranges))
    jj = list(firsts(ranges)) + [len(seq)]
    subs = list(seq[i:j] for i, j in zip(ii, jj))
    return subs


def get_cut_positions(nucs):
    diff_sums = itertools.accumulate(y-x for x, y in nucs)
    cuttings = [0] + list(diff_sums)[:-1]
    return [x-c for (x, _), c in zip(nucs, cuttings)]


def get_args(margin=4):
    parser = argparse.ArgumentParser()
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument(
        "--calls",
        help="gff file with the nucleR calls",
        required=True
    )
    required_named.add_argument(
        "--genome_file",
        help="directory containing the reference genome in fasta format",
        required=True
    )
    required_named.add_argument(
        "--range",
        help="genomic range to process (ex: chrI:2000-2600)",
        required=True
    )
    required_named.add_argument(
        "--seq_output",
        help="output file for the sequence",
        required=True
    )
    required_named.add_argument(
        "--nucs_output",
        help="output file for the nucleosomes",
        required=True
    )
    parser.add_argument(
        "--margin",
        help="number of bases around selected nucleosomes",
        type=int,
        default=margin
    )
    args = parser.parse_args()
    return args.calls, args.genome_file, args.range, args.margin, args.nucs_output, args.seq_output


###############################################################################


def main():
    calls_file, genome_file, range_str, margin, nucs_f, seq_f = get_args()
    chr, start, end = parse_range(range_str)

    ###########################################################################

    with open(calls_file) as fh:
        nucs = list(get_nucs(fh, chr, start, end))

    fasta = read_fasta(genome_file)
    seq = fasta[chr]

    ###########################################################################

    mnucs = merge(nucs)

    sel_start, sel_end = expand_range(get_range(nucs), margin)
    shifted_nucs = [shift(x, sel_start) for x in mnucs]
    sub_seq = seq[sel_start:sel_end]

    nucs_pos = get_cut_positions(shifted_nucs)
    sub_seqs = cut_seq(sub_seq, shifted_nucs)

    ###########################################################################

    out_nucs = ' '.join(map(str, nucs_pos)) + "\n"
    out_seq = ''.join(sub_seqs) + "\n"

    with open(nucs_f, 'w') as fh:
        fh.write(out_nucs)
    with open(seq_f, 'w') as fh:
        fh.write(out_seq)

    return 0

###############################################################################

if __name__ == "__main__":
    sys.exit(main())

###############################################################################
