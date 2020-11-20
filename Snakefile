### Samples ###
ALL_SAMPLES = ["WT", "KO"]

### Define the size of the genomic bins used to generate the HiC contact matrix ###
ALL_BIN_SIZE = [20000, 40000, 150000, 500000, 1000000]

### Define the rules ###
ALL_RAW_MAPS = expand("hic_results/matrix/{sample}/raw/{bsize}/{sample}_{bsize}{suffix}", sample = ALL_SAMPLES, bsize = ALL_BIN_SIZE, suffix = ["_abs.bed", ".matrix"] * len(ALL_BIN_SIZE))
ALL_ICE_MATRIX = expand("hic_results/matrix/{sample}/iced/{bsize}/{sample}_{bsize}_iced.matrix", sample = ALL_SAMPLES, bsize = ALL_BIN_SIZE)


### Rule all ###
rule all:
	input:
		ALL_ICE_MATRIX


### List of rules to be divided into rules .smk format file ###

rule fastp:
	input:
		fastq1="raw/{sample}/{sample}_R1.fastq.gz",
		fastq2="raw/{sample}/{sample}_R2.fastq.gz"
	output:
		trimmed1="fastp/{sample}_R1.fastp.fastq.gz",
		trimmed2="fastp/{sample}_R2.fastp.fastq.gz"
	params:
		parameters="-w 5 -t 1 -A -Q -L"
	log:
		fastpLog="logs/{sample}.fastpLog"
	shell:
		"fastp -i {input.fastq1} -I {input.fastq2} \
		-o {output.trimmed1} -O {output.trimmed2} \
		{params.parameters} \
		2> {log.fastpLog}"


rule alignStep1:
	input:
		aln1="fastp/{sample}_{pos}.fastp.fastq.gz"
	output:
		mapped="mapped_reads/{sample}_{pos}.bam",
		unmapped="unmapped_reads/{sample}_{pos}.unmap.fastq"
	params:
		index="/hpcnfs/techunits/bioinformatics/refdata/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome",
		log="los/{sample}_{pos}.bowtie2.log",
		sm="SM:{sample}_{pos}"
	log:
		alignStep1_log="logs/{sample}_{pos}.alignStep1.log"

	shell:
		"bowtie2 --very-sensitive \
		-L 30 \
		--score-min \
		L,-0.6,-0.2 \
		--end-to-end \
		--reorder \
		--un {output.unmapped} \
		--rg-id BMG --rg \
		{params.sm} -p 3 \
		-x {params.index} \
		-U {input.aln1} \
		2>> {log.alignStep1_log} | \
		samtools view -F 4 -bS - > {output.mapped}"

rule cutsite_trimming:
	input:
		tocut1="unmapped_reads/{sample}_{pos}.unmap.fastq"
	output:
		cut1="unmapped_reads/{sample}_{pos}.unmap.trimmed.fastq"
	params:
		cutSite = "GATC"
	log:
		cutLog="logs/{sample}_{pos}.unmapped_trimming.log"
	shell:
		"/hpcnfs/home/ieo5073/miniconda3/envs/HiCpro/bin/HiC-Pro_2.11.4/scripts/cutsite_trimming \
		--fastq {input.tocut1} \
		--cutsite {params.cutSite} \
		--out {output.cut1} > \
		{log.cutLog} \
		2>&1"

rule alignStep2:
	input:
		aln1="unmapped_reads/{sample}_{pos}.unmap.trimmed.fastq"
	output:
		mapped="mapped_reads/{sample}_{pos}.bwt2glob.unmap_bwt2loc.bam"
	params:
		index="/hpcnfs/techunits/bioinformatics/refdata/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome",
		sm="SM:{sample}_{pos}_genome.bwt2glob.unmap"
	log:
		alignStep2_log="logs/{sample}_{pos}.alignStep2.log"
	shell:
		"bowtie2 --very-sensitive \
		-L 20 \
		--score-min \
		L,-0.6,-0.2 \
		--end-to-end \
		--reorder \
		--rg-id BML --rg \
		{params.sm} -p 3 \
		-x {params.index} \
		-U {input.aln1} \
		2>> {log.alignStep2_log} | \
		samtools view -bS - > {output.mapped}"

rule mappingCombine_merge:
	input:
		bam1="mapped_reads/{sample}_{pos}.bam",
		bam2="mapped_reads/{sample}_{pos}.bwt2glob.unmap_bwt2loc.bam"
	output:
		mergedBam="mapped_reads/bwt2/{sample}/{sample}_{pos}.bwt2merged.bam"
	params:
		thread=6,
		parameters="-n -f"
	shell:
		"/hpcnfs/home/ieo5073/miniconda3/envs/HiCpro/bin/samtools merge \
		-@ {params.thread} \
		{params.parameters} \
		{output.mergedBam} \
		{input.bam1} \
		{input.bam2}"

rule mappingCombine_sort:
	input:
		toSort="mapped_reads/bwt2/{sample}/{sample}_{pos}.bwt2merged.bam"
	output:
		sorted="mapped_reads/bwt2/{sample}/{sample}_{pos}.bwt2merged.sorted.bam"
	params:
		thread=6,
		memory=12,
		parameters="-n -T",
		temporary="tmp/{sample}_{pos}.genome"
	shell:
		"/hpcnfs/home/ieo5073/miniconda3/envs/HiCpro/bin/samtools sort \
		-@ {params.thread} \
		-m {params.memory}G \
		{params.parameters} {params.temporary} \
		-o {output.sorted} \
		{input.toSort}"

rule mergeSAM:
	input:
		toMerge_f="mapped_reads/bwt2/{sample}/{sample}_R1.bwt2merged.sorted.bam",
		toMerge_r="mapped_reads/bwt2/{sample}/{sample}_R2.bwt2merged.sorted.bam"
	output:
		mergedBam="mapped_reads/bwt2/{sample}/{sample}.bwt2pairs.bam"
	params:
		parameters="-q 10 -t -v"
	log:
		mergesam="{sample}.mergeSam.log"
	shell:
		"/hpcnfs/home/ieo5073/miniconda3/envs/HiCpro/bin/python \
		/hpcnfs/home/ieo5073/miniconda3/envs/HiCpro/bin/HiC-Pro_2.11.4/scripts/mergeSAM.py \
		{params.parameters} \
		-f {input.toMerge_f} \
		-r {input.toMerge_r} \
		-o {output.mergedBam} \
		2>{log.mergesam}"

rule mapped_2_hic_fragments:
	input:
		toMap="mapped_reads/bwt2/{sample}/{sample}.bwt2pairs.bam"
	output:
		mapped="hic_results/data/{sample}/{sample}.bwt2pairs.{suffix}"
	params:
		digestedGenome="/hpcnfs/scratch/DP/frossi/ANALYSIS/econway/HiChIP/mm10DpnII_digested",
		parameters="-v -a",
		outdir="hic_results/data/{sample}"
	shell:
		"""
		/hpcnfs/home/ieo5073/miniconda3/envs/HiCpro/bin/python \
		/hpcnfs/home/ieo5073/miniconda3/envs/HiCpro/bin/HiC-Pro_2.11.4/scripts/mapped_2hic_fragments.py \
		{params.parameters} \
		-f {params.digestedGenome} \
		-r {input.toMap} \
		-o {params.outdir}
		
		"""

rule merge_valid_interactions:
	input:
		toMerge="hic_results/data/{sample}/{sample}.bwt2pairs.validPairs"
	output:
		Merged="hic_results/data/{sample}/{sample}.allValidPairs"
	shell:
		"""
		LANG=en; sort -T tmp -S 50% -k2,2V -k3,3n -k5,5V -k6,6n \
		-m {input.toMerge} | \
		awk -F"\t" 'BEGIN{{c1=0;c2=0;s1=0;s2=0}}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){{print;c1=$2;c2=$5;s1=$3;s2=$6}}' > \
		{output.Merged}

		"""

rule build_raw_maps:
	input:
		toBuild="hic_results/data/{sample}/{sample}.allValidPairs"
	output:
		Built="hic_results/matrix/{sample}/raw/{bsize}/{sample}_{bsize}{suffix}"
	params:
		chrSize="/hpcnfs/data/DP/Databases/mm10.chrom.sizes",
		ifile="/dev/stdin",
		m_format="upper",
		bsize_par=ALL_BIN_SIZE,
		Built_par="hic_results/matrix/{sample}/raw/{bsize}/{sample}_{bsize}"
	shell:
		"""
		cat {input.toBuild} | \
		/hpcnfs/home/ieo5073/miniconda3/envs/HiCpro/bin/HiC-Pro_2.11.4/scripts/build_matrix \
		--matrix-format {params.m_format} \
		--binsize {params.bsize_par} \
		--chrsizes {params.chrSize} \
		--ifile {params.ifile} \
		--oprefix {params.Built_par}
		"""

rule ice_normalisation:
	input:
		toNorm="hic_results/matrix/{sample}/raw/{bsize}/{sample}_{bsize}.matrix"
	output:
		Norm="hic_results/matrix/{sample}/iced/{bsize}/{sample}_{bsize}_iced.matrix"
	params:
		filter_Low = 0.02,
		filter_High = 0,
		max_Iter = 100,
		eps = 0.1,
		bias = 1,
		verbose = 1
	log:
		ice_logs = "logs/{sample}_{bsize}_ice.log"
	shell:
		"""
		ice --results_filename {output.Norm} \
		--filter_low_counts_perc {params.filter_Low} \
		--filter_high_counts_perc {params.filter_High} \
		--max_iter {params.max_Iter} \
		--eps {params.eps} \
		--remove-all-zeros-loci \
		--output-bias {params.bias} \
		--verbose {params.verbose} \
		{input.toNorm} \
		>{log.ice_logs}

		"""

	
