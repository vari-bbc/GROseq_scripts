#PBS -l walltime=400:00:00
#PBS -l mem=250gb
#PBS -l nodes=1:ppn=32
#PBS -M ian.beddows@vai.org
#PBS -m abe
#PBS -N GROseq

#for i in {1..17}; do qsub -q bbc ../src/pbs.trimgalore.script.$i; done
#for i in {0..17}; do qsub -q bbc ../src/pbs.star.script.$i; done
#======================================================================= SETUP:
PPN=32 
MEM=250
DIR=/secondary/projects/bbc/research/GROP_20190515_GroSeq/
RAW_DATA=${DIR}raw_data/
TRIMMED_DATA=${DIR}trimmed_data/
#~ STAR_GENOME_DIR=/primary/projects/bbc/references/human/indexes/hg38/hg38_STAR/ # hg38
STAR_GENOME_DIR=/primary/projects/bbc/references/human/indexes/hg38/hg38_STAR_NRSA/ # hg38 with NRSA genome
ref_fasta=${STAR_GENOME_DIR}hg38.fa
ANALYSIS=${DIR}analysis/star/
if [ ! -d ${ANALYSIS} ]; then 
	mkdir ${ANALYSIS}
fi
#======================================================================= Trim data
echo "Generating trimgalore scripts"
files=(${RAW_DATA}*_L000_R1_001.fastq.gz) # 

if [ ! -d ${TRIMMED_DATA} ]; then
	mkdir ${TRIMMED_DATA}
fi

for i in "${!files[@]}"; do
	SAMPLEflowcell1=${files[$i]%%_L000_R1_001.fastq.gz} # Removing the end
	#~ SAMPLEflowcell2=${SAMPLEflowcell1/flowcell1/flowcell2}
	SAMPLE=$(basename $SAMPLEflowcell1)
	
	
cat > ${DIR}src/pbs.trimgalore.script.$i.sh << EOF

#PBS -l walltime=40:00:00
#PBS -l mem=10gb
#PBS -l nodes=1:ppn=1
#PBS -M ian.beddows@vai.org
#PBS -m a
#PBS -N trimgalore_${SAMPLE}

cd ${TRIMMED_DATA}

#~ # Trim for flowcell 1
trim_galore ${SAMPLEflowcell1}_L000_R1_001.fastq.gz \
-q 20 \
--fastqc

# Now trim the trimmed for polyA:
#~ trim_galore ${SAMPLE}_L000_R1_001_trimmed.fq.gz \
#~ -q 20 \
#~ -a "A{51}"
#~ --fastqc

#~ mv ${SAMPLE}_L001_R1_001_trimmed_trimmed.fq.gz ${SAMPLE}_noPolyA_L001_R1_001_trimmed_trimmed.fq.gz
	
EOF

	echo "	${SAMPLE}	pbs.trimgalore.script.$i.sh"

done

#~ exit

#======================================================================= MAPPING:
# Now do the star 2-pass mapping

# use 'find' to find the trimmed files for the samples & give them 
# directly to star
i=0
echo "Generating STAR mapping scripts"
for SAMPLE in `ls ${TRIMMED_DATA}|grep trimmed.fq.gz|sort|uniq`; do # get the unique samples
	SAMPLE=${SAMPLE%_L000_R1_001_trimmed.fq.gz}
	R1_FILES=$(find ${TRIMMED_DATA}|grep _trimmed.fq.gz|grep $SAMPLE\_|tr '\n' ',')
	R1_FILES=$(sed -e 's/\,$//' <<< $R1_FILES) # remove final ,
	
	# For viewing the files:
	# R1_FILES=$(tr ',' '\n' <<< $R1_FILES)
	# echo "$R1_FILES"

cat > ${DIR}src/pbs.star.script.$i.sh << EOF

#PBS -l walltime=40:00:00
#PBS -l mem=250gb
#PBS -l nodes=1:ppn=32
#PBS -M ian.beddows@vai.org
#PBS -m ae
#PBS -N STAR_${SAMPLE}


cd ${ANALYSIS}

#Launch the alignment
STAR --genomeDir ${STAR_GENOME_DIR} \
--readFilesIn \
$R1_FILES \
--twopassMode Basic \
--outReadsUnmapped Fastx \
--runThreadN ${PPN} \
--readFilesCommand zcat \
--outFileNamePrefix ${SAMPLE} \
--quantMode GeneCounts \
--outSAMmapqUnique 60 \
--outSAMtype BAM SortedByCoordinate


# Generate unaligned bam files:
#~ java -Xmx120G -jar /secondary/projects/bbc/tools/picardtools/picard_v2.17.11/picard.jar \
#~ FastqToSam \
#~ FASTQ=$R1_FILES \
#~ SM=${SAMPLE} \
#~ OUTPUT=${SAMPLE}.unaligned.bam


EOF
	echo "$SAMPLE	${DIR}src/pbs.star.script.$i.sh"
	((i++))
done

#~ exit

DEPTHFILE=${DIR}deliverables/depth_of_coverage
if [ ! -f ${DEPTHFILE}.filtered.txt ]; then
	echo "Getting depth data..."
	#~ cat /primary/projects/bbc/references/human/annotation/hg38/ensembl/Homo_sapiens.GRCh38.87.chr.nohead.gtf |cut -f1,3-5|grep gene|cut -f 1,3,4 > ${DIR}hg38.ensembl.gene.bed
	#~ cat /primary/projects/bbc/references/human/annotation/hg38/ensembl/Homo_sapiens.GRCh38.87.chr.nohead.gtf |cut -f3,9|grep -w gene|awk -F ";" '{print $1}'|awk -F "gene_id " '{print $2}'|sed -e 's/\"//g' > ${DIR}hg38.ensembl.geneNames.txt
	samtools depth -q 30 -Q 30 -b ${DIR}hg38.ensembl.gene_plus2kb.bed $(ls -1 ${ANALYSIS}*.bam|grep -v .bai| perl -pe 's/\n/ /g') > ${DEPTHFILE}.txt
	ls -1 ${ANALYSIS}*.bam|grep -v .bai| perl -pe 's/\n/ /g' > ${DIR}header.txt
	#~ paste ${DIR}hg38.ensembl.gene.bed ${DIR}hg38.ensembl.geneNames.txt > ${DIR}hg38.ensembl.gene.txt
	grep -v -P '\d\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0' ${DEPTHFILE}.txt > ${DEPTHFILE}.filtered.txt
fi 


#~ ${DIR}/src/GROP_20190515_GroSeq.pl \
#~ -gtf ${DIR}deliverables/hsapiens_gene_ensembl_plusMinus2kb.txt \
#~ -depth ${DIR}deliverables/depth_of_coverage.filtered.txt \
#~ -outfile ${DIR}Depth_from_TSS_plusMinus.txt

#~ exit


# To get uniq mapping rates to txt file from STAR logs:

if [ ! -f ${DIR}/deliverables/UniquelyMappingRates.txt ]; then
	echo "Generating Mapstats"
	echo "sample	UniquelyMappingRates" > ${DIR}/deliverables/UniquelyMappingRates.txt
	echo "sample	UniquelyMappingReads" > ${DIR}/deliverables/UniquelyMappingReads.txt
	grep 'Uniquely mapped reads %' ${ANALYSIS}*Log.final.out |cut -d ' ' -f1,29|awk -F 'Log.final.out: \\|' '{print $1$2}' >> ${DIR}/deliverables/UniquelyMappingRates.txt
	grep 'Uniquely mapped reads number' ${ANALYSIS}*Log.final.out |cut -d ' ' -f1,24|awk -F 'Log.final.out: \\|' '{print $1$2}' >> ${DIR}/deliverables/UniquelyMappingReads.txt
fi

#~ # to generate unaligned bam files:
if [ ! -f ${DIR}deliverables/joinedSTAR2pass.withheader.count.txt ]; then
	echo "Generating final star matrix:"

	cd ${ANALYSIS}

	# STAR does not output unaligned reads to the file, so no need for -F 4 flag, don't worry about -f/F 256 (multimapped) because only interested in genes & this is RNA data

	# *ReadsPerGene.out.tab:
	# column 1: gene ID
	# column 2: counts for unstranded RNA-seq
	# column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
	# column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

	cat > program.awk << EOF
	BEGIN {FS = "\t"} # field separator
	#~ # {for (i = 2; i <= NF; i += 4) printf ("%s%c", \$i, i + 3 <= NF ? "\t" : "\n");}					# for unstranded
	{for (i = 4; i <= NF; i += 4) printf ("%s%c", \$i, i + 3 <= NF ? "\t" : "\n");}					# for stranded
EOF
	find *ReadsPerGene.out.tab | awk -F 'ReadsPerGene' '{print $1}' |tr '\n' '\t'> header.txt && sed -i -e '$a\' header.txt && sed -i -e 's/^/\t/' header.txt
	paste *ReadsPerGene.out.tab | tail -n +5 | cut -f 1 > rownames
	paste *ReadsPerGene.out.tab | tail -n +5 | awk -f program.awk > joinedSTAR2pass.count.txt
	paste rownames joinedSTAR2pass.count.txt > joinedSTAR2pass.withrownames.count.txt
	cat header.txt joinedSTAR2pass.withrownames.count.txt > joinedSTAR2pass.withheader.count.txt

	ln -s ${ANALYSIS}joinedSTAR2pass.withheader.count.txt ${DIR}deliverables/joinedSTAR2pass.withheader.count.txt
fi
#whether the data is from a strand-specific assay.
#Specify 'yes', 'no', or 'reverse' (default: yes).
#'reverse' means 'yes' with reversed strand interpretation



#======================================================================= MAPPING:
echo "NRSA analysis ..."

# Update R and perl
module load bbc/R/R-3.6.0
export PERL5LIB=$PERL5LIB:/home/ian.beddows/go/GROP_20190515_GroSeq/NRSA-v2/lib/lib/perl5

# First need to convert the bam files to bed format:
# And convert the bam2bed output to 6 column bed format otherwise pause_PROseq.pl does not work


for bam in `ls -1 ${ANALYSIS}*.bam`; do
	
	sample=$(basename $bam);
	sample=${sample%Aligned.sortedByCoord.out.bam}
	echo $sample
	tmp=${ANALYSIS}${sample}.bed.tmp
	bed=${ANALYSIS}${sample}.bed
	if [ ! -f $bed ]; then
		echo "Converting $sample to bed format"
		bam2bed --split < ${bam} > $tmp
		cat $tmp|cut -f1-6 > $bed
		rm $tmp
	else
		echo "bed format for sample $sample exists"
	fi
	
done


cd /secondary/projects/bbc/research/GROP_20190515_GroSeq/analysis/nrsa/


# Now run pause_PROseq for each treatment against the solvent
i=0
declare -a arr=("Com" "PHA" "GFGRO20" "100")
for trt in "${arr[@]}"; do
	ctrl_files=$(find ${ANALYSIS} | grep -E '.bed$'|grep -E 'Sol') #| sed 's/^/-in1 /') # grep Solvent
	trt_files=$(find ${ANALYSIS} | grep -E '.bed$' |grep -E "$trt"|grep -v '\/_')  #| sed 's/^/-in2 /')
	
	outdir=/secondary/projects/bbc/research/GROP_20190515_GroSeq/analysis/nrsa/${trt}_out/
	
	#~ $(echo $gfiles)
	#~ gfiles=$(find ${DE_ASSEMBLY_DIR} | grep -E 'g.vcf$' | sed 's/^/--variant /')
	#~ echo "Files $i"
	#~ echo "	$ctrl_files"
	#~ echo "	$trt_files"
	
	cmd_erna="eRNA.pl -w $outdir -m hg38 -in1 $(echo $ctrl_files) -in2 $(echo $trt_files)"
	cmd="pause_PROseq_DEBUG.pl -o $outdir -m hg38 -in1 $(echo $ctrl_files) -in2 $(echo $trt_files)"
	echo "_________________"
	echo "${trt}"

	if [ ! -d ${outdir}known_gene ]; then
		echo "running the PROseq command"
		echo "${cmd}"
		$cmd
	else
		echo "Directory ${outdir}known_gene for $trt exists .."
	fi
	if [ ! -d ${outdir}eRNA ]; then
		echo "running the eRNA command for $trt"
		echo "${cmd_erna}"
		$cmd_erna
	else
		echo "Directory ${outdir}eRNA for $trt exists .."
	fi
	((i++))
done

#~ exit;

#======================================================================= HOMER MOTIF ANALYSIS GGAA REPEATS:
echo "=========================="
echo "HOMER GGAA REPEAT ANALYSIS"
echo "=========================="
HOMER=/secondary/projects/bbc/tools/homer/bin/
homer=${DIR}analysis/homer/
mkdir -p $homer
cd ${homer}
motif_file=GGAA_REPEAT.motif
if [ ! -f ${motif_file} ]; then
	${HOMER}seq2profile.pl GGAA 0 GGAA_REPEAT > ${motif_file} # create motif with 0 mismatches
fi
motif_result='GGAA_GRCh38.bed'
if [ ! -f $motif_result ]; then
	${HOMER}scanMotifGenomeWide.pl ${motif_file} ${ref_fasta} -bed -mask > ${motif_result} # use the motif file to find these in the ref fasta, masking lowercase
fi
motif_result_merged="merged.${motif_result}"
if [ ! -f $motif_result_merged ]; then
	# now merge the found motifs to within 20bp, count how many intervals get merged with: -c 1 -o count 
	bedtools merge -i ${motif_result} -s -d 20 -c 4,1,6 -o distinct,count,distinct > ${motif_result_merged} # force strandedness, retain chr and strand info
fi
motif_result_merged_filtered=${motif_result_merged/.bed/.filtered.bed}
#echo $motif_result_merged_filtered
if [ ! -f $motif_result_merged_filtered ]; then
	awk '$5>=3' $motif_result_merged > $motif_result_merged_filtered # filter for merged repeats that consist of >=3 individual GGAA repeats
	bedtools sort -chrThenScoreD -i merged.GGAA_GRCh38.filtered.bed > merged.GGAA_GRCh38.filtered.sorted.bed 
fi



#======================================================================= Create stranded bigwig files
ppn=32
source ~/.bashrc
conda activate ucsc_bigwig
bigwig=${DIR}analysis/bigwig/
mkdir -p $bigwig
for bam in `ls ${ANALYSIS}|grep .bam$`; do
	outfile=${bam/Aligned.sortedByCoord.out.bam/.forward.strand.cpm.bw}
	if [ ! -f ${outfile} ]; then 
		echo "generating bigwig for $bam -> ${outfile}"
		bamCoverage -b ${ANALYSIS}${bam} -o ${bigwig}${outfile} -p ${ppn} -of bigwig --normalizeUsing CPM --binSize 10 --filterRNAstrand forward
	fi
	outfile=${bam/Aligned.sortedByCoord.out.bam/.reverse.strand.cpm.bw}
	if [ ! -f ${outfile} ]; then 
		echo "generating bigwig for $bam -> ${outfile}"
		bamCoverage -b ${ANALYSIS}${bam} -o ${bigwig}${outfile} -p ${ppn} -of bigwig --normalizeUsing CPM --binSize 10 --filterRNAstrand reverse
	fi
done







