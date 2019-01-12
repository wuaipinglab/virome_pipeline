#!/usr/bin/env bash

#Make sure you have installed required bioinformatics tools below before running this pipeline
#example parameters (data must be in /data/input/example in compressed fastq format, /data/output directory, /data/database directory and /packge directory and must exist)
#./virome_pipeline.sh /data/input /data/output example type /data/database /package

INPUT_PATH=${1} #path to input directory
OUTPUT_PATH=${2} #path to output directory
SAMPLE_NAME=${3} #Used for directory name of sample
SAMPLE_TYPE=${4} #host type of sample, value = "Culicoides" or "Mosquito"
DATABASE=${5} #Path to reference files and database files directory
PACKAGE=${6} #Path to bioinformatics tools directory

RAW_DATA=${INPUT_PATH}/${SAMPLE_NAME}
SCRIPT=/Pipeline
FASTQC=${PACKAGE}/fastqc
TRIM=${PACKAGE}/Trimmomatic-0.35/trimmomatic-0.35.jar
ADAPTER=${DATABASE}/adapters.fa
STAR=${PACKAGE}/STAR/bin/Linux_x86_64/STAR
MEGAHIT=${PACKAGE}/megahit
BBMAP=${PACKAGE}/bbmap
RSCRIPT=${PACKAGE}/R-3.3.1/bin/Rscript
ORF_FINDER=${PACKAGE}/ORFfinder
CAP3=${PACKAGE}/CAP3/CAP3
BLAST=${PACKAGE}/blast

mkdir ${OUTPUT_PATH}/01qc
mkdir ${OUTPUT_PATH}/02trim
mkdir ${OUTPUT_PATH}/03qc_new
mkdir ${OUTPUT_PATH}/04unmapped_reads
mkdir ${OUTPUT_PATH}/05assembled_contig
mkdir ${OUTPUT_PATH}/06reads_number
mkdir ${OUTPUT_PATH}/07blast_out
mkdir ${OUTPUT_PATH}/08blast_out_FP
mkdir ${OUTPUT_PATH}/09virome
mkdir ${OUTPUT_PATH}/10orf_out
mkdir ${OUTPUT_PATH}/11blastp_out_orf
mkdir ${OUTPUT_PATH}/12orf_annotation

#1.qc and trim
echo "###########################################   Qc and Trim   #####################################################"
#1.1.qc
mkdir ${OUTPUT_PATH}/01qc/${SAMPLE_NAME}
${FASTQC} ${RAW_DATA}/${SAMPLE_NAME}_combined_R1.fastq -t 22 -o ${OUTPUT_PATH}/01qc/${SAMPLE_NAME}
${FASTQC} ${RAW_DATA}/${SAMPLE_NAME}_combined_R2.fastq -t 22 -o ${OUTPUT_PATH}/01qc/${SAMPLE_NAME}
#1.2.trim
mkdir ${OUTPUT_PATH}/02trim/${SAMPLE_NAME}
java -jar ${TRIM} PE -threads 22 ${RAW_DATA}/${SAMPLE_NAME}_combined_R1.fastq ${RAW_DATA}/${SAMPLE_NAME}_combined_R2.fastq -baseout ${OUTPUT_PATH}/02trim/${SAMPLE_NAME}/${SAMPLE_NAME}_combined_trim ILLUMINACLIP:${ADAPTER}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:10
#1.3.qc_again
mkdir ${OUTPUT_PATH}/03qc_new/${SAMPLE_NAME}
${FASTQC} ${OUTPUT_PATH}/02trim/${SAMPLE_NAME}/${SAMPLE_NAME}_combined_trim_1P -t 22 -o ${OUTPUT_PATH}/03qc_new/${SAMPLE_NAME}
${FASTQC} ${OUTPUT_PATH}/02trim/${SAMPLE_NAME}/${SAMPLE_NAME}_combined_trim_2P -t 22 -o ${OUTPUT_PATH}/03qc_new/${SAMPLE_NAME}
echo "${SAMPLE_NAME} qc and trim finished"
echo "#################################################################################################################"

#2.remove host genome and rRNA
echo "#######################################   Remove host genome and rRNA   #########################################"
mkdir ${OUTPUT_PATH}/04unmapped_reads
if [ "${SAMPLE_TYPE}" = "Culicoides" ];then
	#remove Culicoides rRNA
	echo "The host is Culicoides!"
	mkdir ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_rRNA
	$STAR --runMode alignReads --genomeDir ${DATABASE}/Culicoides/rRNA --runThreadN 22 --readFilesIn ${OUTPUT_PATH}/02trim/${SAMPLE_NAME}/${SAMPLE_NAME}_combined_trim_1P ${OUTPUT_PATH}/02trim/${SAMPLE_NAME}/${SAMPLE_NAME}_combined_trim_2P --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMunmapped Within --limitBAMsortRAM 10259511070 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_rRNA/${SAMPLE_NAME}_ --outReadsUnmapped Fastx --outFilterMismatchNmax 60
	#remeove Culicoides genome
	mkdir ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_genome
	$STAR --runMode alignReads --genomeDir ${DATABASE}/Culicoides/genes --runThreadN 22 --readFilesIn ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_rRNA/${SAMPLE_NAME}_Unmapped.out.mate1 ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_rRNA/${SAMPLE_NAME}_Unmapped.out.mate2 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMunmapped Within --limitBAMsortRAM 10259511070 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_genome/${SAMPLE_NAME}_ --outReadsUnmapped Fastx --outFilterMismatchNmax 60
else
	#remove Mosquito rRNA
	echo "The host is Mosquito!"
	mkdir ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_rRNA
	$STAR --runMode alignReads --genomeDir ${DATABASE}/Mosquito/rRNA --runThreadN 22 --readFilesIn ${OUTPUT_PATH}/02trim/${SAMPLE_NAME}/${SAMPLE_NAME}_combined_trim_1P ${OUTPUT_PATH}/02trim/${SAMPLE_NAME}/${SAMPLE_NAME}_combined_trim_2P --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMunmapped Within --limitBAMsortRAM 10259511070 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_rRNA/${SAMPLE_NAME}_ --outReadsUnmapped Fastx --outFilterMismatchNmax 60
	#remove Mosquito genome
	mkdir ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_genome
	$STAR --runMode alignReads --genomeDir ${DATABASE}/Mosquito/genome --runThreadN 22 --readFilesIn ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_rRNA/${SAMPLE_NAME}_Unmapped.out.mate1 ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_rRNA/${SAMPLE_NAME}_Unmapped.out.mate2 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMunmapped Within --limitBAMsortRAM 10259511070 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_genome/${SAMPLE_NAME}_ --outReadsUnmapped Fastx --outFilterMismatchNmax 60
fi
echo "${SAMPLE_NAME} remove host genome and rRNA finished"
echo "#################################################################################################################"

#3.denovo assembly
echo "#######################################   denovo assembly   #####################################################"
${MEGAHIT} -1 ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_genome/${SAMPLE_NAME}_Unmapped.out.mate1 -2 ${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_genome/${SAMPLE_NAME}_Unmapped.out.mate2 -o ${OUTPUT_PATH}/05assembled_contig/${SAMPLE_NAME}_megahit_out --out-prefix ${SAMPLE_NAME}
echo "${SAMPLE_NAME} denovo assembly finished"
echo "#################################################################################################################"

#4.calculate the number of assembled contig
echo "##################################   assembled contig quality evalution   #######################################"
mkdir ${OUTPUT_PATH}/06reads_number/${SAMPLE_NAME}
${BBMAP}/bbwrap.sh ref=${OUTPUT_PATH}/05assembled_contig/${SAMPLE_NAME}_megahit_out/${SAMPLE_NAME}.contigs.fa in=${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_genome/${SAMPLE_NAME}_Unmapped.out.mate1 in2=${OUTPUT_PATH}/04unmapped_reads/${SAMPLE_NAME}_genome/${SAMPLE_NAME}_Unmapped.out.mate2 out=${OUTPUT_PATH}/06reads_number/${SAMPLE_NAME}/${SAMPLE_NAME}.aln.sam.gz local=t maxindel=80 minid=0.8
${BBMAP}/pileup.sh in=${OUTPUT_PATH}/06reads_number/${SAMPLE_NAME}/${SAMPLE_NAME}.aln.sam.gz out=${OUTPUT_PATH}/06reads_number/${SAMPLE_NAME}/${SAMPLE_NAME}_reads_number_out.txt
${RSCRIPT} ${SCRIPT}/read_number.R ${OUTPUT_PATH} ${OUTPUT_PATH}/06reads_number/${SAMPLE_NAME}/${SAMPLE_NAME}_reads_number_out.txt ${OUTPUT_PATH}/06reads_number/${SAMPLE_NAME}/${SAMPLE_NAME}_reads2contig_number.xls
echo "${SAMPLE_NAME} assembled contig quality evalution finished"
echo "#################################################################################################################"

#5.create database and mapped contig with viral nr/nt sequences
echo "###########################################   virus contig   ####################################################"
mkdir ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}
#5.1.make database
#viruses_sequences (vi_nr_seq.fa and vi_nt_seq.fa) were extracted from NCBI nr/nt database according to taxonomy ID
#makeblastdb -in ${DATABASE}/vi_nr/vi_nr_seq.fa -input_type fasta -dbtype prot
#makeblastdb -in ${DATABASE}/vi_nt/vi_nt_seq.fa -input_type fasta -dbtype nucl
#database name:vi_nr_seq.fa vi_nt_seq.fa 
#5.2.mapping
${BLAST}/blastx -query ${OUTPUT_PATH}/05assembled_contig/${SAMPLE_NAME}_megahit_out/${SAMPLE_NAME}.contigs.fa -db ${DATABASE}/vi_nr/vi_nr_seq.fa -evalue 1e-2 -max_target_seqs 1 -num_threads 22 -out ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_out.xls -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen stitle qcovs"
echo "${SAMPLE_NAME} nr ok"
${BLAST}/blastn -query ${OUTPUT_PATH}/05assembled_contig/${SAMPLE_NAME}_megahit_out/${SAMPLE_NAME}.contigs.fa -db ${DATABASE}/vi_nt/vi_nt_seq.fa -evalue 1e-2 -max_target_seqs 1 -num_threads 22 -out ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out.xls -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen stitle qcovs"
echo "${SAMPLE_NAME} nt ok"
${BLAST}/blastx -query ${input_path}/05assembled_contig/${sample_name}_megahit_out/${sample_name}.contigs.fa -db ${refdir}/cdd/cdd_delta -evalue 1e-2 -max_target_seqs 1 -num_threads 22 -out ${input_path}/07blast_out/${sample_name}/${sample_name}_cdd_out.xls -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen stitle qcovs"

sed -i "s/#/ /g" ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_out.xls
perl ${SCRIPT}/cdd_mapped.pl ${DATABASE}/cdd/virus_cdd_UI_all.txt ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_out.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_vi_out.xls
echo "${SAMPLE_NAME} cdd ok"
cat ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_vi_out.xls | awk '{print $1}' | uniq |sort > ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_mapped_contig_id.txt

sed -i "s/#/ /g" ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_out.xls
sed -i "s/#/ /g" ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out.xls
perl ${SCRIPT}/blast_process.pl ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_out.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_out_per.xls
perl ${SCRIPT}/blast_process.pl ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out_per.xls
${RSCRIPT} ${SCRIPT}/nr_nt_length_100_50.R ${OUTPUT_PATH} ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_out_per.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out_per.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_out_100.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out_50.xls

#5.3.merge the nr/nt blast result
cat ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_out_100.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out_50.xls | uniq | sort > ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_contig.xls

#5.4.select viral contig from vi_blast_out
${RSCRIPT} ${SCRIPT}/virus_nr_mapped_contig.R ${OUTPUT_PATH} ${DATABASE}/taxaNodes_taxaNames.RData ${OUTPUT_PATH}/06reads_number/${SAMPLE_NAME}/${SAMPLE_NAME}_reads2contig_number.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_contig.xls ${DATABASE}/accessionTaxa.sql ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_mapped_contig_taxonomy_information.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_mapped_virus_contig_taxonomy_information.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_mapped_virus_contig_id.txt
perl ${SCRIPT}/find_mapped_seq.pl ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_mapped_virus_contig_id.txt ${OUTPUT_PATH}/05assembled_contig/${SAMPLE_NAME}_megahit_out/${SAMPLE_NAME}.contigs.fa ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_mapped_unique_contigs.fa
echo "${SAMPLE_NAME} vi_mapped contig ok"
echo "#################################################################################################################"

#6.remove false-positive contig
echo "###########################################   remove false-positive contig   ####################################"
#6.1.mapping against nr db
mkdir ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}
${BLAST}/blastn -query ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_mapped_unique_contigs.fa -db ${DATABASE}/nt/nt -evalue 1e-10 -max_target_seqs 1 -num_threads 22 -out ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out_FP.xls -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen stitle qcovs"
sed -i "s/#/ /g" ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out_FP.xls
mkdir ${OUTPUT_PATH}/09virome/${SAMPLE_NAME}
#6.2.remove fp contig
${RSCRIPT} ${SCRIPT}/remove_false_positive_contig.R ${OUTPUT_PATH} ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out_FP.xls ${DATABASE}/taxaNodes_taxaNames.RData ${DATABASE}/accessionTaxa.sql ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_mapped_contig_taxonomy_information.xls ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_false_positive_contig_taxonomy_information.xls ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_false_positive_contig_id.txt ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_mapped_virus_contig_taxonomy_information.xls ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_remove_fp_contig_taxonomy_information.xls ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_remove_fp_contig_id.txt
${RSCRIPT} ${SCRIPT}/virome.R ${OUTPUT_PATH} ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_remove_fp_contig_taxonomy_information.xls ${SAMPLE_NAME} ${OUTPUT_PATH}/09virome/${SAMPLE_NAME}/${SAMPLE_NAME}_family.xls ${OUTPUT_PATH}/09virome/${SAMPLE_NAME}/${SAMPLE_NAME}_virome_family.pdf ${OUTPUT_PATH}/09virome/${SAMPLE_NAME}/${SAMPLE_NAME}_species.xls ${OUTPUT_PATH}/09virome/${SAMPLE_NAME}/${SAMPLE_NAME}_virome_species.pdf
#6.3.family_seq
perl ${SCRIPT}/find_mapped_seq.pl ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_remove_fp_contig_id.txt ${OUTPUT_PATH}/05assembled_contig/${SAMPLE_NAME}_megahit_out/${SAMPLE_NAME}.contigs.fa ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_remove_fp_contig_seq.fa
${RSCRIPT} ${SCRIPT}/family_seq.R ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME} ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_remove_fp_contig_taxonomy_information.xls
for i in $(ls ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/*_family_id.txt)
do
	family=`basename $i |sed 's/_family_id.txt//'`
	perl ${SCRIPT}/find_mapped_seq.pl ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_id.txt ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_remove_fp_contig_seq.fa ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_seq.fa
	${cap3}/cap3 ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_seq.fa > ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family.log
	perl ${SCRIPT}/family_id.pl ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_seq.fa.cap.contigs ${family} ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_seq.fa.cap.contigs.new
	perl ${SCRIPT}/family_id.pl ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_seq.fa ${family} ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_seq.fa.new
	cat ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_seq.fa.cap.contigs.new ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_seq.fa.new > ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${family}_family_seq.cap.final.fa	
done
cat ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/*_family_seq.cap.final.fa > ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/all_family_seq.cap.final.fa
echo "${SAMPLE_NAME} remove false-positive contig ok"
echo "#################################################################################################################"

#7.ORFfinder
echo "###########################################   ORF prediction   ##################################################"
mkdir ${OUTPUT_PATH}/10orf_out/${SAMPLE_NAME}
perl ${SCRIPT}/find_mapped_seq.pl ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_mapped_contig_id.txt ${OUTPUT_PATH}/05assembled_contig/${SAMPLE_NAME}_megahit_out/${SAMPLE_NAME}.contigs.fa ${OUTPUT_PATH}/10orf_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_contig_seq.fa
cat ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/all_family_seq.cap.final.fa ${OUTPUT_PATH}/10orf_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_contig_seq.fa > ${OUTPUT_PATH}/10orf_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_nt_all_family_contig_seq.fa
${ORF_FINDER} -in ${OUTPUT_PATH}/10orf_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_nt_all_family_contig_seq.fa -out ${OUTPUT_PATH}/10orf_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_nt_all_family_orf.fa
echo "${SAMPLE_NAME} orf prediction ok"
echo "#################################################################################################################"

#8.ORF annotation by blastp against nr database
echo "###########################################   ORF blastp   ######################################################"
mkdir ${OUTPUT_PATH}/11blastp_out_orf/${SAMPLE_NAME}
${BLAST}/blastp -query ${OUTPUT_PATH}/10orf_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_nt_all_family_orf.fa -db ${refdir}/nr/nr -evalue 1e-5 -max_target_seqs 1 -num_threads 22 -out ${OUTPUT_PATH}/11blastp_out_orf/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_nt_all_family_orf_blast.xls -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen stitle qcovs"
sed -i "s/#/ /g" ${OUTPUT_PATH}/11blastp_out_orf/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_nt_all_family_orf_blast.xls
echo "${SAMPLE_NAME} orf blastp ok"
echo "#################################################################################################################"

#9.collect ORF blastp information,including:ORF number,ORF length,ORF sequence,blast_out
echo "###########################################   ORF annotation   ##################################################"
mkdir ${OUTPUT_PATH}/12orf_annotation/${SAMPLE_NAME}
${RSCRIPT} ${SCRIPT}/orf_annotation.R ${OUTPUT_PATH} ${OUTPUT_PATH}/11blastp_out_orf/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_nt_all_family_orf_blast.xls ${DATABASE}/accessionTaxa.sql ${DATABASE}/taxaNodes_taxaNames.RData ${OUTPUT_PATH}/12orf_annotation/${SAMPLE_NAME}/${SAMPLE_NAME}_orf_annotation.xls ${OUTPUT_PATH}/12orf_annotation/${SAMPLE_NAME}/${SAMPLE_NAME}_virus_related_orf_annotation.xls
echo "${SAMPLE_NAME} orf annotation ok"												
echo "#################################################################################################################"

#10.merge the vi.nr/vi.nt/nt/ cdd/ orf/ blast result
echo "###########################################   merge result   ####################################################"
mkdir ${OUTPUT_PATH}/13merge_out/${SAMPLE_NAME}
cat ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_out_per.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out_per.xls ${OUTPUT_PATH}/08blast_out_FP/${SAMPLE_NAME}/${SAMPLE_NAME}_nt_out_FP.xls | sort > ${OUTPUT_PATH}/13merge_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_nt_blast_out.xls
${RSCRIPT} ${SCRIPT}/taxonomizr.R ${OUTPUT_PATH}/13merge_out/${SAMPLE_NAME} ${DATABASE}/taxaNodes_taxaNames.RData ${DATABASE}/accessionTaxa.sql ${OUTPUT_PATH}/13merge_out/${SAMPLE_NAME}/${SAMPLE_NAME}_nr_nt_blast_out.xls ${OUTPUT_PATH}/07blast_out/${SAMPLE_NAME}/${SAMPLE_NAME}_cdd_out.xls ${OUTPUT_PATH}/12orf_annotation/${SAMPLE_NAME}/${SAMPLE_NAME}_orf_annotation.xls ${SAMPLE_NAME}_nr_nt_blast_out_taxonomy.xls ${SAMPLE_NAME}_cdd_out.xls ${SAMPLE_NAME}_orf_annotation.xls ${OUTPUT_PATH}/13merge_out/${SAMPLE_NAME}/${SAMPLE_NAME}_merge_nr_nt_cdd_orf_blast_out_taxonomy.xls
echo "#################################################################################################################"




