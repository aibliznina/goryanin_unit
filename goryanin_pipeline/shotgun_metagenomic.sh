# Shotgun metagenomic analysis pipeline including 6 step of analysis processes
# 0) Checking qualitiy of seqeuncing data
# 1) Filtering good quality reads of each metatrascriptomic
# 2) Assigning taxonomic of each metagenomic data
# 3) Assembling read to contig and scaffold
# 4) Calling gene and extract protein seqeunce from assembled read
# 5) Annoting protein sequence
# 6) Comparing fuctional of each metagenomic

# Load configuration file
source configuration

# Load fastq file
list_fa_file=();
# echo "Load fastq file"
# for entry in `find $dir_seq -type f -name *.fastq*`; do
#     list_fa_file+=($entry)
# done
list_fa_file=(
	"GO-electricigeus-1_S1_L001_R1_001.fastq,GO-electricigeus-1_S1_L001_R2_001.fastq" 
	"GO-electricigeus-2_S2_L001_R1_001.fastq,GO-electricigeus-2_S2_L001_R2_001.fastq" 
	"GO-electricigeus-3_S3_L001_R1_001.fastq,GO-electricigeus-3_S3_L001_R2_001.fastq"
)
# list_fa_short_name=(
# 	"GO01S01L001"
# 	"GO02S02L001"
# 	"GO03S03L001"
# )

# STEP 0: Checking quality of seqeuncing data
dir_out_qc=0_seq_qc
mkdir -p dir_out_qc
for entry in ${list_fa_file[@]}; do
	fastq_paired_end=(${entry//,/ })
	echo "${dir_fastqc} ${dir_seq}${fastq_paired_end[0]} -t ${cpu_max} &"
	${dir_fastqc} ${dir_seq}${fastq_paired_end[0]} -t ${cpu_max} &
	echo "${dir_fastqc} ${dir_seq}${fastq_paired_end[1]} -t ${cpu_max}"
	${dir_fastqc} ${dir_seq}${fastq_paired_end[1]} -t ${cpu_max}
done
mv ${dir_seq}*fastqc.* dir_out_qc

# STEP 1: Filtering good quality reads of each metatrascriptomic
# trimmed left and right read that low quality (>= 20 bp) and maximum N is one base on each read.
# output in directory seq_filtered and sufix file is ${entry_name}.1_ns
dir_out_fread=1_seq_filtered
mkdir -p ${dir_out_fread}
for entry in ${list_fa_file[@]}; do
	fastq_paired_end=(${entry//,/ })
	entry_name=(${entry//_R1_/ })
	entry_name=${entry_name[0]}
	echo perl ${dir_printseq} -verbose -fastq ${dir_seq}${fastq_paired_end[0]} -fastq2 ${dir_seq}${fastq_paired_end[1]} -ns_max_n 0 -out_good ${dir_out_fread}/${entry_name}.1_ns -out_bad ${dir_out_fread}/${entry_name}.with_ns -trim_qual_left 20 -trim_qual_right 20
	perl ${dir_printseq} -verbose -fastq ${dir_seq}${fastq_paired_end[0]} -fastq2 ${dir_seq}${fastq_paired_end[1]} -ns_max_n 0 -out_good ${dir_out_fread}/${entry_name}.1_ns -out_bad ${dir_out_fread}/${entry_name}.with_ns -trim_qual_left 20 -trim_qual_right 20 &
done

# STEP 2: Assigning taxonomy of each metagenomic data
# Output of this step is taxonomy tab-separation files of each metagenomic data
dir_out_taxa=2_taxonomy
mkdir -p ${dir_out_taxa}
mkdir -p ${dir_out_taxa}/visualize_taxonomy
for entry in ${list_fa_file[@]}; do
	fastq_paired_end=(${entry//,/ })
	entry_name=(${entry//_R1_/ })
 	entry_name=${entry_name[0]}
	
	echo ${dir_metaxa} -1 ${dir_out_fread}/${entry_name}.1_ns_1.fastq -2 ${dir_out_fread}/${entry_name}.1_ns_2.fastq --quality_filter=F --quality_trim=T --save_raw=T -o ${dir_out_taxa}/${entry_name} --cpu=${cpu_max}
	${dir_metaxa} -1 ${dir_out_fread}/${entry_name}.1_ns_1.fastq -2 ${dir_out_fread}/${entry_name}.1_ns_2.fastq --quality_filter=F --quality_trim=T --save_raw=T -o ${dir_out_taxa}/${entry_name} --cpu=${cpu_max}
	
	mkdir -p ${dir_out_taxa}/${entry_name}
	mv ${dir_out_taxa}/${entry_name}.* ${dir_out_taxa}/${entry_name}/
	mkdir -p ${dir_out_taxa}/${entry_name}/fasta_marker_16s
	mv ${dir_out_taxa}/${entry_name}/*.fasta ${dir_out_taxa}/${entry_name}/fasta_marker_16s/

	# count the number of species, genera by Metaxa2 Taxonomic Traversal Tool (metaxa2_ttt)
	${dir_metaxa2_ttt} -i ${dir_out_taxa}/${entry_name}/$entry_name.taxonomy.txt -o ${dir_out_taxa}/${entry_name}/${entry_name}
	mkdir -p ${dir_out_taxa}/${entry_name}/taxonomy
	mv ${dir_out_taxa}/${entry_name}/*.level* ${dir_out_taxa}/${entry_name}/taxonomy/

	# taxonomy visualization
	cp ${dir_out_taxa}/${entry_name}/taxonomy/$entry_name.level_${taxonomy_level}.txt ${dir_out_taxa}/visualize_taxonomy
	awk -F \\t '{print $2 FS $1}' ${dir_out_taxa}/visualize_taxonomy/${entry_name}.level_${taxonomy_level}.txt | sed 's/;/\t/g' > ${dir_out_taxa}/visualize_taxonomy/${entry_name}.level_${taxonomy_level}.converted.txt
	ktImportText ${dir_out_taxa}/visualize_taxonomy/${entry_name}.level_${taxonomy_level}.converted.txt -o ${dir_out_taxa}/visualize_taxonomy/${entry_name}.level_${taxonomy_level}.html
done
# LSU rRNA
dir_out_taxa=2_taxonomy_LSU
mkdir -p ${dir_out_taxa}
mkdir -p ${dir_out_taxa}/visualize_taxonomy
for entry in ${list_fa_file[@]}; do
	fastq_paired_end=(${entry//,/ })
	entry_name=(${entry//_R1_/ })
 	entry_name=${entry_name[0]}
	
	echo ${dir_metaxa} -1 ${dir_out_fread}/${entry_name}.1_ns_1.fastq -2 ${dir_out_fread}/${entry_name}.1_ns_2.fastq --quality_filter=F --quality_trim=T --save_raw=T -o ${dir_out_taxa}/${entry_name} --cpu=${cpu_max} -g lsu
	${dir_metaxa} -1 ${dir_out_fread}/${entry_name}.1_ns_1.fastq -2 ${dir_out_fread}/${entry_name}.1_ns_2.fastq --quality_filter=F --quality_trim=T --save_raw=T -o ${dir_out_taxa}/${entry_name} --cpu=${cpu_max} -g lsu
	
	mkdir -p ${dir_out_taxa}/${entry_name}
	mv ${dir_out_taxa}/${entry_name}.* ${dir_out_taxa}/${entry_name}/
	mkdir -p ${dir_out_taxa}/${entry_name}/fasta_marker_16s
	mv ${dir_out_taxa}/${entry_name}/*.fasta ${dir_out_taxa}/${entry_name}/fasta_marker_16s/

	# count the number of species, genera by Metaxa2 Taxonomic Traversal Tool (metaxa2_ttt)
	${dir_metaxa2_ttt} -i ${dir_out_taxa}/${entry_name}/$entry_name.taxonomy.txt -o ${dir_out_taxa}/${entry_name}/${entry_name}
	mkdir -p ${dir_out_taxa}/${entry_name}/taxonomy
	mv ${dir_out_taxa}/${entry_name}/*.level* ${dir_out_taxa}/${entry_name}/taxonomy/

	# taxonomy visualization
	cp ${dir_out_taxa}/${entry_name}/taxonomy/$entry_name.level_${taxonomy_level}.txt ${dir_out_taxa}/visualize_taxonomy
	awk -F \\t '{print $2 FS $1}' ${dir_out_taxa}/visualize_taxonomy/${entry_name}.level_${taxonomy_level}.txt | sed 's/;/\t/g' > ${dir_out_taxa}/visualize_taxonomy/${entry_name}.level_${taxonomy_level}.converted.txt
	ktImportText ${dir_out_taxa}/visualize_taxonomy/${entry_name}.level_${taxonomy_level}.converted.txt -o ${dir_out_taxa}/visualize_taxonomy/${entry_name}.level_${taxonomy_level}.html
done

# STEP3: Assembling read to contig and scaffold
dir_out_assem=3_assembled
mkdir -p dir_out_assem
for entry in ${list_fa_file[@]}; do
	fastq_paired_end=(${entry//,/ })
	entry_name=(${entry//_R1_/ })
	echo ${dir_SPAdes} --pe1-1 1_seq_filtered/${entry_name}.1_ns_1.fastq --pe1-2 1_seq_filtered/${entry_name}.1_ns_2.fastq -o ${dir_out_assem}/${entry_name}-assembled --meta --threads ${cpu_max}
	${dir_SPAdes} --pe1-1 1_seq_filtered/${entry_name}.1_ns_1.fastq --pe1-2 1_seq_filtered/${entry_name}.1_ns_2.fastq -o ${dir_out_assem}/${entry_name}-assembled --meta --threads ${cpu_max}
done

# # STEP4: Calling gene and extract protein seqeunce from assembled read
dir_out_cgene=4_gene_call
mkdir -p ${dir_out_cgene}
for entry in ${list_fa_file[@]}; do
	fastq_paired_end=(${entry//,/ })
	entry_name=(${entry//_R1_/ })
	echo ${dir_FragGeneScan} -s dir_out_assem/${entry_name}.fa -o ${dir_out_cgene}/${entry_name}-calling_gene -w 0 -t illumina_5 -p ${cpu_max}
	${dir_FragGeneScan} -s dir_out_assem/${entry_name}.fa -o ${dir_out_cgene}/${entry_name}-calling_gene -w 0 -t illumina_5 -p ${cpu_max}

	# Change fasta head name
	echo python change_gene_name.py -i ${dir_out_cgene}/${entry_name}-calling_gene.faa -a ${dir_out_cgene}/${entry_name}-calling_gene.out -o ${dir_out_cgene}/${entry_name}.prot.faa
	python change_gene_name.py -i ${dir_out_cgene}/${entry_name}-calling_gene.faa -a ${dir_out_cgene}/${entry_name}-calling_gene.out -o ${dir_out_cgene}/${entry_name}.prot.faa &

done

# STEP5: Annoting protein sequence divied into 2 part
dir_out_anno=5_annotation
mkdir -p ${dir_out_anno}
# PART A Annotating of each protein by GhostKOALA (http://www.kegg.jp/ghostkoala/) by manually
#        Input file is protein sequence from gene calling tool in directory "${dir_out_cgene}" in ${entry_name}-calling_gene-contigs.faa
#        Output download from websites

# PART B Annotating protein domain by InterProScan
mkdir -p ${dir_out_anno}/temp_fasta/
mkdir -p ${dir_out_anno}/protein_domain_prediction/
mkdir -p ${dir_out_anno}/log_interproscan

for entry in ${list_fa_file[@]}; do
	fastq_paired_end=(${entry//,/ })
	entry_name=(${entry//_R1_/ })

	# Split fasta file
	echo python ${dir_split_fasta} -i ${dir_out_cgene}/${entry_name}.prot.faa -o ${dir_out_anno}/temp_fasta --node 30
	python ${dir_split_fasta} -i ${dir_out_cgene}/${entry_name}.prot.faa -o ${dir_out_anno}/temp_fasta --node 30
	
	# # sbatch runing on each node
	work_dir=${dir_out_anno}/temp_fasta/${entry_name}/
	list_sh_file=();
	for entry in `find $work_dir -type f -name run_interpro*`; do
	    list_sh_file+=($entry)
	done
	for sh in ${list_sh_file[@]}; do
	    echo sbatch $sh
	    sbatch $sh
	done
	
	echo mv ${entry_name}.prot.faa.* ${dir_out_anno}/protein_domain_prediction/
	mv mv ${entry_name}.prot.faa.* ${dir_out_anno}/protein_domain_prediction/
	echo mv *_InterProScan_* ${dir_out_anno}/log_interproscan/
	mv mv *_InterProScan_* ${dir_out_anno}/log_interproscan/

done


# # STEP6: Functional comparative between microbial communities

# # STEP_additional: bining contig, scaffold, and read
# # Convert a paired-end read FastQ file to fasta file
# mkdir -p dir_out_assem_bining
# for entry in ${list_fa_file[@]}; do
# 	fastq_paired_end=(${entry//,/ })
# 	entry_name=(${entry//_R1_/ })

# 	echo python ${dir_fq2fa} -1 ${dir_out_fread}/${entry_name}.1_ns_1.fastq -2 ${dir_out_fread}/${entry_name}.1_ns_2.fastq -o dir_out_assem_bining/${entry_name}.fa
# 	# python ${dir_fq2fa} -1 ${dir_out_fread}/${entry_name}.1_ns_1.fastq -2 ${dir_out_fread}/${entry_name}.1_ns_2.fastq -o dir_out_assem_bining/${entry_name}.fa

# 	echo ${dir_MetaAnnotator} dir_out_assem_bining/${entry_name}.fa dir_out_assem/${entry_name}-assembled/contigs.fasta ${dir_MetaAnnotator_bin}outfilename1.idx ${dir_MetaAnnotator_bin}nodes.dmp ${dir_MetaAnnotator_bin}outfilename1.idx.mtx ${dir_MetaAnnotator_bin}outfilename2.idx
# 	${dir_MetaAnnotator} dir_out_assem_bining/${entry_name}.fa dir_out_assem/${entry_name}-assembled/contigs.fasta ${dir_MetaAnnotator_bin}outfilename1.idx ${dir_MetaAnnotator_bin}nodes.dmp ${dir_MetaAnnotator_bin}outfilename1.idx.mtx ${dir_MetaAnnotator_bin}outfilename2.idx
# done

