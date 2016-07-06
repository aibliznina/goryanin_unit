import sys, os, stat
from optparse import OptionParser

# Parse Command Line
usage = """
Description: tool for convert name of each fasta file
convert header of each protein sequence name to systematic name
for example 
>Manes.S111100.1.p pacid=32327512 transcript=Manes.S111100.1 locus=Manes.S111100 ID=Manes.S111100.1.v6.1 annot-version=v6.1
1)Organism name
2)Running number
3)transcript name
4)locus name
5)ID name
6)annotation version
"""

parser = OptionParser(usage=usage)
parser.add_option("-i", "--input", dest="dir_input",
	default=None, metavar="FILE",
	help="FASTA file")
parser.add_option("-o", "--output", dest="dir_out",
	default=None, metavar="FILE",
	help="Output fasta file name")
parser.add_option("-s", "--sample", dest="sample",
	default=None, type="string",
	help="Organism name or sample name for show in header of each fasta" )
options,args = parser.parse_args()

if not options.dir_input:
	sys.exit("Missing input input fasta file, -i input.fa")
if not os.path.isfile(options.dir_input):
	sys.exit("Missing input input fasta file %r" % options.dir_input)
if not options.dir_out:
	sys.exit("Missing output fasta file name, -o output.fa")
dir_input = os.path.abspath(options.dir_input)
dir_out_fa = options.dir_out
sample = options.sample

fasta_input = open(dir_input).read().split('>')
fasta_input = fasta_input[1:]

w_file_fa_out = open(dir_out_fa, 'w')
w_file_postion = open(dir_out_fa + ".pos.txt", 'w')
w_file_postion.write("contig_name\tstrand\tstart\tend\n")

catalog_contig = {}

num = 0
for i in fasta_input[:40]:
	header = i[:i.find('\n')].split('_')
	contig_name = header[0] + "_" + header[1]
	if contig_name not in catalog_contig.keys() :
		catalog_contig[contig_name] = 1
	else:
		catalog_contig[contig_name] += 1
	strand = header[-1]
	start = header[-3]
	end = header[-2]
	sequence = i[i.find('\n'):-1].replace('\n', '')
	seq_len = len(sequence)
	w_file_postion.write(contig_name + "\t" + strand + "\t" + start + "\t" + end + "\n")
	seq_name = ">" + sample + "."+contig_name +"P%06d" % catalog_contig[contig_name] + ".p" 
	print(seq_name, seq_len)
	w_file_fa_out.write(seq_name + "\n" + sequence + "\n")
	num += 1

# for key,value in catalog_contig:
# 	print( key, value)

