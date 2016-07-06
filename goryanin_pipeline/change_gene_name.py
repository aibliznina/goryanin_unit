import sys, os, stat
from optparse import OptionParser

# Parse Command Line
usage = """
Description: tool for convert gene name result from FragGeneScan
convert header of each protein sequence name to systematic name
for example 
>GO1S1L001.con.000000_62_370_+
>GO1S1L001.con.000000.000000
"""

parser = OptionParser(usage=usage)
parser.add_option("-i", "--input_fasta", dest="dir_input_fasta",
	default=None, metavar="FILE",
	help="FASTA file")
parser.add_option("-a", "--input_annotation", dest="dir_input_ann",
	default=None, metavar="FILE",
	help="tab-delimeted file")
parser.add_option("-o", "--output", dest="dir_out",
	default=None, metavar="FILE",
	help="Output fasta file name")
options,args = parser.parse_args()

if not options.dir_input_fasta:
	sys.exit("Missing input input fasta file, -i input.fa")
if not options.dir_input_ann:
	sys.exit("Missing input input annotation file, -a ann.out")
if not os.path.isfile(options.dir_input_fasta):
	sys.exit("Missing input input fasta file %r" % options.dir_input_fasta)
if not os.path.isfile(options.dir_input_ann):
	sys.exit("Missing input input annotation file %r" % options.dir_input_ann)
if not options.dir_out:
	sys.exit("Missing output fasta file name, -o output.fa")
dir_input_fasta = os.path.abspath(options.dir_input_fasta)
dir_input_ann = os.path.abspath(options.dir_input_ann)
dir_out_fa = options.dir_out

fasta_input = open(dir_input_fasta).read().split('>')
fasta_input = fasta_input[1:]

annotation_input = open(dir_input_ann).read().split('>')
annotation_input = annotation_input[1:]

w_file_fa_out = open(dir_out_fa, 'w')
w_file_fa_log = open(dir_out_fa + '.log', 'w')

catalog_protein = {}

num = 0
for i in fasta_input:
	header = i[:i.find('\n')]
	seq = i[i.find('\n'):].replace('\n','')
	catalog_protein[header] = seq

for i in annotation_input:
	header = i[:i.find('\n')]
	info = i[i.find('\n')+1:-1].split('\n')
	num_protein_in_contig = 0
	if(info[0] != ''):
		for j in info:
			j = j.split('\t')
			fasta_header = header + '_' + j[0] + '_' + j[1] + '_' + j[2]
			if fasta_header in catalog_protein.keys():
				new_header = header + '.' + "%06d" % num_protein_in_contig
				# print('>' + new_header)
				w_file_fa_out.write('>' + new_header + '\n' + catalog_protein[fasta_header] + '\n')
				w_file_fa_log.write(new_header + '\t' + fasta_header + '\n')
				num_protein_in_contig += 1
			else:
				print("Error not found this sequence in fasta file.")
				exit()

