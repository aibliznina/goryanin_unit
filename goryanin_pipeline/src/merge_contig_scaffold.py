import sys, os, stat
from optparse import OptionParser

def getGCcontent(sequence):
	GC = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
	AT = sequence.count('A') + sequence.count('T') + sequence.count('a') + sequence.count('t')
	return float(GC) * 100 / (AT + GC)

# Parse Command Line
usage = """
Description: tool for convert name of each fasta file
convert header of each protein sequence name to systematic name
for example 
>GO1_1.Cont.000001
>GO1_1.Scaf.000001
>sample.Cont/Scaf
1) Sample
2) Contig or scaffold
3) running number
"""

parser = OptionParser(usage=usage)
parser.add_option("-c", "--contig", dest="dir_contig",
	default=None, metavar="FILE",
	help="FASTA file")
parser.add_option("-s", "--scaffold", dest="dir_scaffold",
	default=None, metavar="FILE",
	help="FASTA file")
parser.add_option("-o", "--output", dest="dir_out",
	default=None, metavar="FILE",
	help="Output fasta file name")
parser.add_option("-a", "--sample", dest="sample",
	default=None, type="string",
	help="Organism name or sample name for show in header of each fasta")
options, args = parser.parse_args()

if not options.dir_contig:
	sys.exit("Missing input input contig fasta file, -c contig.fa")
if not os.path.isfile(options.dir_contig):
	sys.exit("Missing input input contig fasta file %r" % options.dir_contig)
if not options.dir_out:
	sys.exit("Missing output fasta file name, -o output.fa")
dir_contig = os.path.abspath(options.dir_contig)
dir_out_fa = options.dir_out
sample = options.sample

fasta_contig = open(dir_contig).read().split('>')
fasta_contig = fasta_contig[1:]

w_file_fa_out = open(dir_out_fa, 'w')
w_file_postion = open(dir_out_fa + ".pos.txt", 'w')
w_file_postion.write("contig_name\told_name\tlength\tGC_content\tnumber_of_N\n")

catalog_contig = {}

num = 0
for i in fasta_contig:
	header = i[:i.find('\n')]
	contig_name = sample + ".con." + "%06d" % num
	sequence = i[i.find('\n'):-1].replace('\n', '')
	seq_len = len(sequence)

	# print(contig_name, seq_len)
	w_file_postion.write(contig_name + "\t" + header + "\t" + str(seq_len) + "\t" + "%.2f" % getGCcontent(sequence) + "\t" + str(sequence.count('N')+ sequence.count('n')) + "\n") 
	w_file_fa_out.write(">"+ contig_name + "\n" + sequence + "\n") 
	num += 1

if options.dir_scaffold:
	fasta_scaffold = open(options.dir_scaffold).read().split('>')
	fasta_scaffold = fasta_scaffold[1:]
	num = 0
	for i in fasta_scaffold:
		header = i[:i.find('\n')]
		contig_name = sample + ".sca." + "%06d" % num
		sequence = i[i.find('\n'):-1].replace('\n', '')
		seq_len = len(sequence)

		# print(contig_name, seq_len)
		w_file_postion.write(contig_name + "\t" + header + "\t" + str(seq_len) + "\t" + "%.2f" % getGCcontent(sequence) + "\t" + str(sequence.count('N') + sequence.count('n')) + "\n") 
		w_file_fa_out.write(">"+ contig_name + "\n" + sequence + "\n")
		num += 1