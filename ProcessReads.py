#Ryan A. Melnyk
#schmelnyk@gmail.com
#UBC Microbiology - Haney Lab

import os
import argparse
from Bio import SeqIO
import subprocess

def parse_args():
	parser = argparse.ArgumentParser(description='''
Pipeline for parsing STAMP sequencing data. Provide path to the sequencing directory as a single argument.
	''')
	parser.add_argument('seqdata', type=str,help='path to the sequencing data')
	parser.add_argument('outdir',type=str,help='path to create output folder')
	return parser.parse_args()

def unzip_files(seqdata,outdir):
	for f in os.listdir(seqdata):
		filename = os.listdir(os.path.join(seqdata,f))[0]
		cmd = 'gunzip {}'.format(os.path.join(seqdata,f,filename))
		os.system(cmd)
		filename = os.listdir(os.path.join(seqdata,f))[0]
		cmd = 'cp {} {}'.format(os.path.join(seqdata,f,filename),os.path.join(outdir,"raw",filename.split("_")[0]+'.fastq'))
		os.system(cmd)
	return

def filter_by_length(f,outdir):
	fq = open(os.path.join(outdir,"filtered",f.split(".")[0]+".filtered.fq"),'w')
	for s in SeqIO.parse(open(os.path.join(outdir,"raw",f),'r'),'fastq'):
		if len(s.seq) == 62:
			SeqIO.write(s,fq,'fastq')
	return

def trim_tn(f,outdir):
	cmds = "scythe -a conserved_region.fna {}.filtered.fq -o {}.trimmed.fq -q sanger".format(os.path.join(outdir,"filtered",f.split(".")[0]),os.path.join(outdir,"trimmed",f.split(".")[0]))
	proc = subprocess.Popen(cmds.split())
	proc.wait()
	return

def get_lengths(f,outdir):
	o = open(os.path.join(outdir,"exact",f.split(".")[0]+".exact.fq"),'w')
	count = 0
	total = 0
	for s in SeqIO.parse(open(os.path.join(outdir,"trimmed",f.split(".")[0]+".trimmed.fq"),'r'),'fastq'):
		total += 1
		if len(s.seq) == 30:
			count += 1
			SeqIO.write(s,o,'fastq')
	return

def find_unique_barcodes(f,outdir):
	barcodes = {}
	for s in SeqIO.parse(open(os.path.join(outdir,"exact",f.split(".")[0]+".exact.fq"),'r'),'fastq'):
		if str(s.seq) not in barcodes:
			barcodes[str(s.seq)] = 1
		else:
			barcodes[str(s.seq)] += 1
	o = open(os.path.join(outdir,"barcodes",f.split(".")[0]+".barcodes.txt"),'w')
	for b in sorted(barcodes,key = lambda x: barcodes[x],reverse=True):
		o.write("{}\t{}\n".format(b, barcodes[b]))
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
	for f in ["raw","filtered","trimmed","exact","barcodes"]:
		try:
			os.makedirs(os.path.join(outdir,f))
		except OSError:
			pass
	unzip_files(seqdata,outdir)

	for f in os.listdir(os.path.join(outdir,"raw")):
		filter_by_length(f,outdir)
		trim_tn(f,outdir)
		get_lengths(f,outdir)
		find_unique_barcodes(f,outdir)

if __name__ == '__main__':
	main()
