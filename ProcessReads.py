# Ryan A. Melnyk
# schmelnyk@gmail.com

import os
import sys
import argparse
from Bio import SeqIO
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(description='''
Pipeline for parsing STAMP sequencing data. Provide path to the sequencing
directory as a single argument.
    ''')
    parser.add_argument('seqdata', type=str,
                        help='path to the sequencing data')
    parser.add_argument('outdir', type=str,
                        help='path to create output folder')
    return parser.parse_args()


def unzip_files(seqdata, outdir):
    for f in os.listdir(seqdata):
        filename = os.listdir(os.path.join(seqdata, f))[0]
        cmd = 'gunzip {}'.format(os.path.join(seqdata, f, filename))
        os.system(cmd)
        filename = os.listdir(os.path.join(seqdata, f))[0]
        cmd = 'cp {} {}'.format(os.path.join(seqdata, f, filename),
                                os.path.join(outdir,
                                             "raw",
                                             filename.split("_")[0]+'.fastq'))
        os.system(cmd)
    return


def filter_by_length(f, outdir):
    fq = open(os.path.join(outdir,
                           "filtered",
                           f.split(".")[0]+".filtered.fq"), 'w')
    for s in SeqIO.parse(open(os.path.join(outdir, "raw", f), 'r'), 'fastq'):
        if len(s.seq) == 62:
            SeqIO.write(s, fq, 'fastq')
    return


def trim_tn(f, outdir):
    cmds = ["scythe",
            "-a",
            os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])),
                         "conserved_region.fna"),
            "-M", "5",
            "{}.filtered.fq".format(os.path.join(outdir,
                                                 "filtered",
                                                 f.split(".")[0])),
            "-o",
            "{}.trimmed.fq".format(os.path.join(outdir,
                                                "trimmed",
                                                f.split(".")[0])),
            "-q",
            "sanger"]
    proc = subprocess.Popen(cmds)
    proc.wait()
    return


def get_lengths(f, outdir):
    o = open(os.path.join(outdir, "exact", f.split(".")[0]+".exact.fq"), 'w')
    count = 0
    total = 0
    for s in SeqIO.parse(open(
            os.path.join(outdir, "trimmed", f.split(".")[0]+".trimmed.fq"),
            'r'),
            'fastq'):
        total += 1
        if len(s.seq) == 30:
            count += 1
            SeqIO.write(s, o, 'fastq')
    return


def find_unique_barcodes(f, outdir):
    barcodes = {}
    for s in SeqIO.parse(open(os.path.join(outdir, "exact",
                         f.split(".")[0]+".exact.fq"), 'r'), 'fastq'):
        if str(s.seq) not in barcodes:
            barcodes[str(s.seq)] = 1
        else:
            barcodes[str(s.seq)] += 1
    o = open(os.path.join(outdir,
                          "barcodes",
                          f.split(".")[0]+".barcodes.txt"), 'w')
    for b in sorted(barcodes, key=lambda x: barcodes[x], reverse=True):
        o.write("{}\t{}\n".format(b, barcodes[b]))
    o.close()
    return


def get_stats(outdir):
    stats = {}
    o = open(os.path.join(outdir, "read_stats.txt"), 'w')
    o.write("sample\traw_reads\texact\n")
    for f in os.listdir(os.path.join(outdir, "raw")):
        count = 0
        for seq in SeqIO.parse(open(os.path.join(outdir, "raw", f),
                               'r'), 'fastq'):
            count += 1
        stats[f.split(".")[0]] = {"raw": count}
        exact_count = 0
        for seq in SeqIO.parse(open(os.path.join(outdir,
                                                 "exact",
                                                 f.split(".")[0] +
                                                 ".exact.fq"), 'r'), 'fastq'):
            exact_count += 1
        stats[f.split(".")[0]]["exact"] = exact_count
        o.write("{}\t{}\t{}\n".format(f.split(".")[0],
                                      str(count),
                                      str(exact_count)))
    o.close()
    return


def make_unique_fasta(f, outdir):
    o = open(os.path.join(outdir, "unique", f.split(".")[0]+".unique.fa"), 'w')
    barcodes = {}
    for line in open(os.path.join(outdir,
                                  "barcodes",
                                  f.split(".")[0]+".barcodes.txt"), 'r'):
        seq = line.rstrip().split("\t")[0]
        o.write(">{}\n".format(seq))
        o.write("{}\n".format(seq))
        barcodes[seq] = int(line.rstrip().split("\t")[1])
    o.close()
    return barcodes


def uclust(f, outdir):
    cmds = "vsearch --cluster_fast {} --minseqlength 30 --uc {} --id 0.8"

    proc = subprocess.Popen(
        cmds.format(os.path.join(outdir,
                                 "unique",
                                 f.split(".")[0]+".unique.fa"),
                    os.path.join(outdir,
                                 "uc",
                                 f.split(".")[0]+".uc")).split())
    proc.wait()
    return


def parse_uclust(f, outdir, barcodes):
    clusters = {}
    for line in open(os.path.join(outdir, "uc", f.split(".")[0]+".uc"), 'r'):
        vals = line.rstrip().split()
        if vals[0] == "#":
            continue
        elif vals[0] == "S":
            clusters[vals[8]] = []
        elif vals[0] == "H":
            clusters[vals[9]].append(vals[8])
        else:
            pass
    o = open(os.path.join(outdir, "clusters", f.split(".")[0]+".clusters.txt"),
             'w')
    for c in clusters:
        match_count = 0
        for match in clusters[c]:
            match_count += barcodes[match]
        o.write("{}\t{}\t{}\n".format(c, barcodes[c], match_count))
    o.close()
    return clusters


def cluster_matrix(outdir):
    all_data = {}
    names = []
    for f in os.listdir(os.path.join(outdir, "clusters")):
        name = f.split(".")[0]
        names.append(name)
        for line in open(os.path.join(outdir, "clusters", f), 'r'):
            vals = line.rstrip().split("\t")
            if vals[0] not in all_data:
                all_data[vals[0]] = {name: int(vals[1]) + int(vals[2])}
            else:
                all_data[vals[0]][name] = int(vals[1]) + int(vals[2])

    o = open(os.path.join(outdir, "barcode_matrix.txt"), 'w')
    o.write("barcode\t{}\n".format("\t".join(names)))
    for x in all_data:
        this_line = []
        for n in names:
            if n in all_data[x]:
                this_line.append(all_data[x][n])
            else:
                this_line.append(0)
        o.write("{}\t{}\n".format(x, "\t".join([str(z) for z in this_line])))
    o.close()
    return


def main():
    args = parse_args()
    seqdata = os.path.abspath(args.seqdata)
    outdir = os.path.abspath(args.outdir)
    try:
        os.makedirs(outdir)
    except OSError:
        pass
    for f in ["raw", "filtered", "trimmed", "exact", "barcodes", "unique",
              "uc", "clusters"]:
        try:
            os.makedirs(os.path.join(outdir, f))
        except OSError:
            pass
    unzip_files(seqdata, outdir)

    for f in os.listdir(os.path.join(outdir, "raw")):
        filter_by_length(f, outdir)
        trim_tn(f, outdir)
        get_lengths(f, outdir)
        find_unique_barcodes(f, outdir)
        barcodes = make_unique_fasta(f, outdir)
        uclust(f, outdir)
        clusters = parse_uclust(f, outdir, barcodes)

    cluster_matrix(outdir)
    get_stats(outdir)


if __name__ == '__main__':
    main()
