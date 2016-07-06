import sys, os, stat
from optparse import OptionParser

# Parse Command Line
usage = """
Description: tool for split fasta file

"""

parser = OptionParser(usage=usage)
parser.add_option("-i", "--input", dest="dir_input",
	default=None, metavar="FILE",
	help="protein fasta file")
parser.add_option("-o", "--output", dest="dir_out",
	default=None, metavar="FILE",
	help="Output filename (required)")
parser.add_option("-n", "--node", dest="node",
	default=None, type="int",
	help="Output filename (required)")

options,args = parser.parse_args()

if not options.dir_input:
	sys.exit("Missing input fasta file, -i xxx.faa")
if not os.path.isfile(options.dir_input):
	sys.exit("Missing input fasta file, %r" % options.dir_input)
if not options.dir_out:
	sys.exit("Missing output directory, -o xxxxx")
if not options.node:
	sys.exit("Missing number of node")

dir_input = os.path.abspath(options.dir_input)
dir_pre_input, input_file_name = os.path.split(dir_input)
fasta_name = input_file_name.split('.')[0]

dir_out = os.path.abspath(options.dir_out)

node = options.node

# Create directory collection of splited fasta files
dir_working=dir_out+"/"+fasta_name
if not os.path.exists(dir_working):
	os.makedirs(dir_working)

fasta = open(dir_input).read().replace('*','').split('>')
fasta = fasta[1:]

numberOfProtSeq = len(fasta)
# print ("Number of sequence: " + str(numberOfProtSeq))
fasta_split_per_file = int(numberOfProtSeq/node) + 1
# print (fasta_split_per_file)
i = 0
file_num = 0
while i < numberOfProtSeq:
	if i % fasta_split_per_file == 0:
		w_file_name = dir_working + "/" + input_file_name + ".%04d" % file_num
		# print ("Write file: " + w_file_name)
		w_fasta = open(w_file_name, 'w')

		# Generate file bash for run InterProScan
		file_name_sh = dir_working + "/run_interpro_" + "%04d" % file_num + ".sh"
		# print("Write file: " + file_name_sh)
		w_sh_run_interpro = open(file_name_sh, 'w')
		
		txt = "#!/bin/bash" + \
			"\n#SBATCH --job-name=" + "%02d" % file_num + "Ips_" + fasta_name + \
			"\n#SBATCH --partition=compute" + \
			"\n#SBATCH --time=24:00:00" + \
			"\n#SBATCH --mem=120g" + \
			"\n#SBATCH -c 24" + \
			"\n#SBATCH --ntasks=1" + \
			"\n#SBATCH --mail-user=nattawet.sriwichai@oist.jp" + \
			"\n#SBATCH --mail-type=END" + \
			"\n#SBATCH --input=none" + \
			"\n#SBATCH --output=job_%j_InterProScan_" + "%04d" % file_num +".out" + \
			"\n#SBATCH --error=job_%j_InterProScan_" + "%04d" % file_num +".err" + \
			"\n\nmodule load java-jdk/1.8.0_20" + \
			"\n\ndir_InterProScan=$HOME/bin/software/annotation/interproscan-5.18-57.0/interproscan.sh" + \
			"\ncpu_max=24" + \
			"\n\necho $dir_InterProScan -i " + w_file_name + " -t p -dp -pa --goterms --iprlookup" + \
			"\n$dir_InterProScan -i " + w_file_name + " -t p -dp -pa --goterms --iprlookup\n" 
		w_sh_run_interpro.write(txt)
		w_sh_run_interpro.close()

		# Change mode to excuteable file
		st = os.stat(file_name_sh)
		os.chmod(file_name_sh, st.st_mode | stat.S_IEXEC)

		file_num += 1
	header = fasta[i][:fasta[i].find('\n')]
	sequence = fasta[i][fasta[i].find('\n'):-1].replace('\n', '')
	seq_len = len(sequence)
	w_fasta.write(">"+header+"\n"+sequence+"\n")
	# print(header, "\t", seq_len)
	i += 1

