#!/bin/bash
#SBATCH --job-name=assembly
#SBATCH --partition=compute
#SBATCH --time=64:00:00
#SBATCH --mem=125g
#SBATCH -c 24
#SBATCH --ntasks=1
#SBATCH --mail-user=nattawet.sriwichai@oist.jp
#SBATCH --mail-type=END
#SBATCH --input=none
#SBATCH --output=job_%j_assembly.out
#SBATCH --error=job_%j_assembly.err

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
dir_out_qc=${dir_out}/${dir_out_qc}
dir_out_fread=${dir_out}/${dir_out_fread}
dir_out_taxa=${dir_out}/${dir_out_taxa}
dir_out_assem=${dir_out}/${dir_out_assem}
dir_out_bining=${dir_out}/${dir_out_bining}
dir_out_cgene=${dir_out}/${dir_out_cgene}
dir_out_anno=${dir_out}/${dir_out_anno}



# # Load fastq file
list_fa_file=()
list_fa_short_name=()
echo "Please type short name"
for R1 in `cd ${dir_seq} && ls *_R1_*`; do
	fastq_paired_end=(${R1} ${R1//_R1_/_R2_})
	case "$R1" in
		*.gz ) 
			sample_name=${R1%%.*}
			;;
		* )
			sample_name=${R1%.*}
			;;
	esac
	sample_name=${sample_name//_R1_001/}
	sample_name=${sample_name//_/}
	sample_name=${sample_name//-/}
	# read -e -p " + ${R1}: " -i $sample_name sample_name
	list_fa_short_name+=(${sample_name})
	list_fa_file+=(${fastq_paired_end[0]},${fastq_paired_end[1]})
done
echo "--------------------------------------------"

# list_fa_file=("GO-EL-01_S1_L001_R1_001.fastq,GO-EL-01_S1_L001_R2_001.fastq" "GO-EL-02_S2_L001_R1_001.fastq.gz,GO-EL-02_S2_L001_R2_001.fastq.gz")
# list_fa_file=("GO-EL-02_S2_L001_R1_001.fastq.gz,GO-EL-02_S2_L001_R2_001.fastq.gz")
# list_fa_file=("${list_fa_file[@]:1}")

echo " ___________________________________________________________"
echo "|         FastQ file name          |      Short name        |"
echo "|__________________________________|________________________|"
num_list=${#list_fa_file[@]}
for (( i=0; i<${#list_fa_file[@]}; i++ )); do
	echo "  ${list_fa_file[$i]%,*}     ${list_fa_short_name[$i]}  "
done
echo "|__________________________________|________________________|"

echo "Shotgun metagenomics analyis pipeline including: "
if [ "${step_0_check_quality}" == "Y" ]; then
	echo "+ Run step 0: Check sequence quality"
fi
if [ "${step_1_filtering}" == "Y" ]; then
	echo "+ Run step 1: Filtering read"
fi
if [ "${step_2_taxonomy_assigning}" == "Y" ]; then
	echo "+ Run step 2: Taxonomic assigenment"
fi
if [ "${step_3_assembling}" == "Y" ]; then
	echo "+ Run step 3: Assembly metagneomic"
fi
if [ "${step_3_1_coassembling}" == "Y" ]; then
	echo "+ Run step 3.1: coassembly metagenomic sample"
fi
if [ "${step_3_2_reassembling}" == "Y" ]; then
	echo "+ Run step 3.2: reassembly read"
fi
if [ "${step_4_gene_calling}" == "Y" ]; then
	echo "+ Run step 4: Gene calling"
fi
if [ "${step_5_annotating}" == "Y" ]; then
	echo "+ Run step 5: gene annotation"
fi




# # STEP 0: Checking quality of seqeuncing data
if [ "${step_0_check_quality}" == "Y" ]; then
	echo "========== STEP 0: Checking quality of seqeuncing data =========="
	mkdir -p ${dir_out_qc}
	for entry in ${list_fa_file[@]}; do
		fastq_paired_end=(${entry//,/ })
		echo "${dir_fastqc} ${dir_seq}/${fastq_paired_end[0]} --threads ${cpu_max} --outdir ${dir_out_qc}"
		${dir_fastqc} ${dir_seq}/${fastq_paired_end[0]} --threads ${cpu_max} --outdir ${dir_out_qc} &
		echo "${dir_fastqc} ${dir_seq}/${fastq_paired_end[1]} --threads ${cpu_max} --outdir ${dir_out_qc}"
		${dir_fastqc} ${dir_seq}/${fastq_paired_end[1]} --threads ${cpu_max} --outdir ${dir_out_qc}
	done
	echo "================================================================"
	rm ${dir_out_qc}/*.zip
fi

# # STEP 1: Filtering good quality reads of each metatrascriptomic
# # trimmed left and right read that low quality (>= 20 bp) and maximum N is one base on each read.
# # output in directory seq_filtered and sufix file is sample name
if [ "${step_1_filtering}" == "Y" ]; then
	echo "========== STEP 1: Filtering good quality reads of each metatrascriptomic =========="
	mkdir -p ${dir_out_fread}
	mkdir -p ${dir_out_fread}/log
	mkdir -p ${dir_out_fread}/qc
	num=0
	for entry in ${list_fa_file[@]}; do
		fastq_paired_end=(${entry//,/ })

		echo "Filtering: ${list_fa_short_name[num]}"
		# Trimmomatic filtering read
		java -jar ${dir_trimmomatic} PE \
			-threads 24 \
			-trimlog ${dir_out_fread}/log/${list_fa_short_name[num]}.log \
			${dir_seq}/${fastq_paired_end[0]} \
			${dir_seq}/${fastq_paired_end[1]} \
			${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R1.unpaired.fq.gz \
			${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R2.unpaired.fq.gz \
			ILLUMINACLIP:/home/n/nattawet-sriwichai/bin/software/filter/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
			LEADING:${p_trimming_LEADING} \
			TRAILING:${p_trimming_TRAILING} \
			SLIDINGWINDOW:${p_trimming_SLIDINGWINDOW} \
			MINLEN:${p_trimming_MINLEN} > ${dir_out_fread}/log/${list_fa_short_name[num]}.summary.txt 2>&1

		# Merge filtered read
		echo "------------------- MERGE -------------------"
		cat ${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R1.unpaired.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq.gz
		cat ${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R2.unpaired.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}_R2.fq.gz
		# rm ${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R1.unpaired.fq.gz
		# rm ${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R2.unpaired.fq.gz
		# Quality checking
		cat ${dir_out_fread}/${list_fa_short_name[num]}_R1.unpaired.fq.gz ${dir_out_fread}/${list_fa_short_name[num]}_R2.unpaired.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}.unpaired.fq.gz
		echo "-------------- Recheck quality --------------"
		echo "Check quality of filtered read: ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq.gz"
		${dir_fastqc} ${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz -t ${cpu_max} --outdir ${dir_out_fread}/qc &
		${dir_fastqc} ${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz -t ${cpu_max} --outdir ${dir_out_fread}/qc
		
		echo "---------------Check number of read --------------"
		num_read_R1=$(cat ${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz | echo $((`wc -l`/4)))
		num_read_R2=$(cat ${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz | echo $((`wc -l`/4)))
		num_read_unpaired=$(cat ${dir_out_fread}/${list_fa_short_name[num]}.unpaired.fq.gz | echo $((`wc -l`/4)))
		printf "${list_fa_short_name[num]}\t${num_read_R1}\t${num_read_R2}\t${num_read_unpaired}\n" >> num_read.txt
		
		num=$((num+1))
	done
	echo "================================================================"
fi


# # STEP 2: Assigning taxonomy of each metagenomic data
# # Output of this step is taxonomy tab-separation files of each metagenomic data
if [ "${step_2_taxonomy_assigning}" == "Y" ]; then
	echo "========== STEP 2: Assigning taxonomy of each metagenomic data =========="
	mkdir -p ${dir_out_taxa}
	mkdir -p ${dir_out_taxa}/visualize_taxonomy
	dir_out_taxa2=${dir_out_taxa}/MetaPlnAn
	mkdir -p ${dir_out_taxa2}
	mkdir -p ${dir_out_taxa2}/visualize_taxonomy

	num=0
	for entry in ${list_fa_file[@]}; do
		fastq_paired_end=(${entry//,/ })
		gunzip -c ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq
		gunzip -c ${dir_out_fread}/${list_fa_short_name[num]}_R2.fq.gz > ${dir_out_fread}/${list_fa_short_name[num]}_R2.fq
		
		${dir_metaxa} \
			-1 ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq \
			-2 ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq \
			-o ${dir_out_taxa}/${list_fa_short_name[num]} \
			-t ${p_taxonomy_organism} \
			-T ${p_taxonomy_identity} \
			-m ${p_taxonomy_match} \
			--cpu=${cpu_max}

		rm ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq ${dir_out_fread}/${list_fa_short_name[num]}_R2.fq
		
		mkdir -p ${dir_out_taxa}/${list_fa_short_name[num]}
		mv ${dir_out_taxa}/${list_fa_short_name[num]}.* ${dir_out_taxa}/${list_fa_short_name[num]}/
		mkdir -p ${dir_out_taxa}/${list_fa_short_name[num]}/fasta_marker_16s
		mv ${dir_out_taxa}/${list_fa_short_name[num]}/*.fasta ${dir_out_taxa}/${list_fa_short_name[num]}/fasta_marker_16s/

		# count the number of species, genera by Metaxa2 Taxonomic Traversal Tool (metaxa2_ttt)
		${dir_metaxa2_ttt} -i ${dir_out_taxa}/${list_fa_short_name[num]}/${list_fa_short_name[num]}.taxonomy.txt \
			-o ${dir_out_taxa}/${list_fa_short_name[num]}/${list_fa_short_name[num]}
		mkdir -p ${dir_out_taxa}/${list_fa_short_name[num]}/taxonomy
		mv ${dir_out_taxa}/${list_fa_short_name[num]}/*.level* ${dir_out_taxa}/${list_fa_short_name[num]}/taxonomy/

		# taxonomy visualization
		cp ${dir_out_taxa}/${list_fa_short_name[num]}/taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.txt ${dir_out_taxa}/visualize_taxonomy
		awk -F \\t '{print $2 FS $1}' ${dir_out_taxa}/visualize_taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.txt | sed 's/;/\t/g' > ${dir_out_taxa}/visualize_taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.converted.txt
		ktImportText ${dir_out_taxa}/visualize_taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.converted.txt \
			-o ${dir_out_taxa}/visualize_taxonomy/${list_fa_short_name[num]}.level_${p_taxonomy_level}.html
		
		# # STEP2.1 high confident of Taxonomy identification based on specific marker gene by MetaPlnAn2: 
		echo "Run STEP2.1 high confident of Taxonomy identification based on specific marker gene by MetaPlnAn2: ${list_fa_short_name[num]}"
	 	${dir_MetaPlnAn} ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq.gz,${dir_out_fread}/${list_fa_short_name[num]}_R2.fq.gz \
	 		--bowtie2out ${dir_out_taxa2}/${list_fa_short_name[num]}.bowtie2.bz2 \
	 		--nproc ${cpu_max} \
	 		--input_type fastq \
	 		--ignore_viruses \
	 		--ignore_eukaryotes \
	 		--min_cu_len 2000 \
	 		-t rel_ab_w_read_stats \
	 		> ${dir_out_taxa2}/${list_fa_short_name[num]}.profile.txt
		num=$((num+1))
	done
	echo "================================================================"
fi

# # STEP3: Assembling read to contig and scaffold
if [ "${step_3_assembling}" == "Y" ]; then
	mkdir -p ${dir_out_assem}
	contig_file_name=""
	num=0
	for entry in ${list_fa_file[@]}; do
		echo "Assembling: ${list_fa_short_name[num]}"
		fastq_paired_end=(${entry//,/ })
		
		${dir_SPAdes} \
			--meta \
			--threads ${cpu_max} \
			--pe1-1 ${dir_out_fread}/${list_fa_short_name[num]}_R1.paired.fq.gz \
			--pe1-2 ${dir_out_fread}/${list_fa_short_name[num]}_R2.paired.fq.gz \
			--pe1-s ${dir_out_fread}/${list_fa_short_name[num]}.unpaired.fq.gz \
			-o ${dir_out_assem}/${list_fa_short_name[num]}
		
		# Change heading of each contig and/or merge contig with scaffold
		python src/merge_contig_scaffold.py \
			--contig ${dir_out_assem}/${list_fa_short_name[num]}/contigs.fasta \
			--output ${dir_out_assem}/${list_fa_short_name[num]}.fa \
			--sample ${list_fa_short_name[num]} 

		contig_file_name="${contig_file_name}${dir_out_assem}/${list_fa_short_name[num]}.fa "
		num=$((num+1))
	done

	# # STEP3.1: Check quality of metagenomic assembly by using metaQuast
	echo "Check quality of metagenomic assembly"
	${dir_metaQuast} ${contig_file_name} -o ${dir_out_assem}/quality --threads=${cpu_max} --max-ref-number 1000

fi

# # STEP 3.2 Co-assembly between data
if [ "${step_3_1_coassembling}" == "Y" ]; then
	mkdir -p ${dir_out_co_assem}
	if [ ${#list_fa_file[@]} -gt 9 ]; then
		echo "metagenomic data set more than 9 data"
		num_data=9
	else
		num_data=${#list_fa_file[@]}
	fi

	text=""
	for (( i=0; i<${num_data}; i++ )); do
		echo "Assembling: ${list_fa_short_name[i]}"
		text="${text}--pe$((i+1))-1 ${dir_out_fread}/${list_fa_short_name[i]}_R1.paired.fq.gz \\
			--pe$((i+1))-2 ${dir_out_fread}/${list_fa_short_name[i]}_R2.paired.fq.gz \\
			--pe$((i+1))-s ${dir_out_fread}/${list_fa_short_name[i]}.unpaired.fq.gz \\
			"
	done

	echo -e "${dir_SPAdes} \\
		--meta \\
		--threads ${cpu_max} \\
		${text} \\
		-o ${dir_out_co_assem}"

	# # # STEP3.1: Check quality of metagenomic assembly by using metaQuast
	# echo "Check quality of metagenomic assembly"
	# ${dir_metaQuast} ${contig_file_name} -o ${dir_out_assem}/quality --threads=${cpu_max} --max-ref-number 1000

fi


# # STEP 3.3: Bining and reassembly
if [ "${step_3_2_reassembling}" == "Y" ]; then
	dir_out_bining=3_bining
	mkdir ${dir_out_bining}

	num=1
	for entry in ${list_fa_file[@]}; do
		fastq_paired_end=(${entry//,/ })

		${dir_MaxBin} -contig "${contig_file_name}${dir_out_assem}/${list_fa_short_name[num]}.fa" \
			-reads ${dir_out_fread}/${list_fa_short_name[num]}_R1.fq.gz \
			-reads2 ${dir_out_fread}/${list_fa_short_name[num]}_R2.fq.gz \
			-out ${dir_out_bining}/${list_fa_short_name[num]} \
			-thread ${cpu_max} \
			-reassembly
		num=$((num+1))
	done
fi


# # STEP4: Calling gene and extract protein seqeunce from assembled read
if [ "${step_4_gene_calling}" == "Y" ]; then
	mkdir -p ${dir_out_cgene}
	num=0
	for entry in ${list_fa_file[@]}; do
		fastq_paired_end=(${entry//,/ })

		echo ${dir_FragGeneScan} -s ${dir_out_assem}/${list_fa_short_name[num]}.fa -o ${dir_out_cgene}/${list_fa_short_name[num]}-calling_gene -w 0 -t illumina_5 -p ${cpu_max}
		${dir_FragGeneScan} -s ${dir_out_assem}/${list_fa_short_name[num]}.fa \
			-o ${dir_out_cgene}/${list_fa_short_name[num]}-calling_gene \
			-w 0 \
			-t illumina_5 \
			-p ${cpu_max}

		# Change fasta head name
		python src/change_gene_name.py -i ${dir_out_cgene}/${list_fa_short_name[num]}-calling_gene.faa \
			-a ${dir_out_cgene}/${list_fa_short_name[num]}-calling_gene.out \
			-o ${dir_out_cgene}/${list_fa_short_name[num]}.prot.faa
		num=$(($num+1))
	done
fi

# # STEP5: Annoting protein sequence divied into 2 part
if [ "${step_5_annotating}" == "Y" ]; then
	mkdir -p ${dir_out_anno}
	# PART A Annotating of each protein by GhostKOALA (http://www.kegg.jp/ghostkoala/) by manually
	#        Input file is protein sequence from gene calling tool in directory "${dir_out_cgene}" in ${entry_name}-calling_gene-contigs.faa
	#        Output download from websites

	# PART B Annotating protein domain by InterProScan
	mkdir -p ${dir_out_anno}/temp_fasta/
	mkdir -p ${dir_out_anno}/protein_domain_prediction/
	mkdir -p ${dir_out_anno}/log_interproscan

	num=0
	for entry in ${list_fa_file[@]}; do

		fastq_paired_end=(${entry//,/ })


		# Split fasta file
		echo python ${dir_split_fasta} -i ${dir_out_cgene}/${list_fa_short_name[num]}.prot.faa -o ${dir_out_anno}/temp_fasta --node 30
		python ${dir_split_fasta} -i ${dir_out_cgene}/${list_fa_short_name[num]}.prot.faa -o ${dir_out_anno}/temp_fasta --node 30
		
		# # sbatch runing on each node
		work_dir=${dir_out_anno}/temp_fasta/${list_fa_short_name[num]}/
		list_sh_file=();
		for entry in `find $work_dir -type f -name run_interpro*`; do
		    list_sh_file+=($entry)
		done
		for sh in ${list_sh_file[@]}; do
		    echo sbatch $sh
		    sbatch $sh
		done
		
		echo mv ${list_fa_short_name[num]}.prot.faa.* ${dir_out_anno}/protein_domain_prediction/
		mv ${list_fa_short_name[num]}.prot.faa.* ${dir_out_anno}/protein_domain_prediction/
		echo mv *_InterProScan_* ${dir_out_anno}/log_interproscan/
		mv *_InterProScan_* ${dir_out_anno}/log_interproscan/
		num=$(($num+1))
	done
fi

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

