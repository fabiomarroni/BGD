#!/bin/bash

#Variables you need to change to your local setting
###################################################
FUNCDIR=/full/path/to/functions/
BASEDIR=/base/directory/
###################################################

#Download database (no need to uncompress it, dada2 can handle gzipped DBs)
#This is needed only once
cd ${BASEDIR}/databases/silva
#Taxonomy
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
#Species assignment. Not sure we will use it.
wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz




#######################################################
#
#  ______ _                  _                _            _                 
# |  ____| |                | |              | |          (_)                
# | |__  | | __ ___   _____ | |__   __ _  ___| |_ ___ _ __ _ _   _ _ __ ___  
# |  __| | |/ _` \ \ / / _ \| '_ \ / _` |/ __| __/ _ \ '__| | | | | '_ ` _ \ 
# | |    | | (_| |\ V / (_) | |_) | (_| | (__| ||  __/ |  | | |_| | | | | | |
# |_|    |_|\__,_| \_/ \___/|_.__/ \__,_|\___|\__\___|_|  |_|\__,_|_| |_| |_|
#                                                                           
########################################################


#Extract all Flavobacterium species from kraken-generated reads
#1. Get all flavobacterium species we classified
cat ${BASEDIR}/05_kraken2/*_S.bracken.txt | cut -f1,2 | grep "Flavobacterium" | sort | uniq | awk -F "\t" '{print "kraken:taxid|"$2}' > ${BASEDIR}/01_data/Flavob_list.txt
#2. Loop over generated reads and extract flavobacterium
READDIR=/projects/marroni/BGD/05_kraken2
OUTDIR=${BASEDIR}/06_onlyflavoreads
for READ1 in ${READDIR}/*.fastq.gz
do
#echo $READ1
OUTFILE=${READ1/$READDIR/$OUTDIR}
OUTFILE=${OUTFILE/.fastq.gz/}
echo Extracting reads to $OUTFILE ...
zcat $READ1 | grep -f ${BASEDIR}/01_data/Flavob_list.txt -A 3 | grep -v '^--'  > ${OUTFILE}.fastq 
done
#3. zip files
gzip ${OUTDIR}/*fastq

#Align with bwa-mem on the Flavobacterium glacei reference from IZSV
#######################################################
#  ______          __          __  __ ______ __  __  
# |  _ \ \        / /\        |  \/  |  ____|  \/  | 
# | |_) \ \  /\  / /  \ ______| \  / | |__  | \  / | 
# |  _ < \ \/  \/ / /\ \______| |\/| |  __| | |\/| | 
# | |_) | \  /\  / ____ \     | |  | | |____| |  | | 
# |____/   \/  \/_/    \_\    |_|  |_|______|_|  |_| 
#                                                    
#######################################################
                                                    

#FULL REFSEQ+BISELLI on all reads (we provide the custom Ref_plus_biselli.fasta file)
bwa-mem2 index ${BASEDIR}/02_reference/Ref_plus_biselli.fasta
#We need the new samtools 

READDIR=${BASEDIR}/00_reads
ALDIR=$BASEDIR/08_bwamem_all_RefSeq_biselli
NTHREAD=2
for READ1 in ${READDIR}/*_R1_001.fastq.gz
do
#echo $READ1
READ2=${READ1/_R1_001.fastq.gz/_R2_001.fastq.gz}
FNAME=$(basename $READ1)
OUTFILE=${READ1/$READDIR/$ALDIR}
OUTFILE=${OUTFILE/_L001_R1_001.fastq.gz/}
if [ -s ${OUTFILE}.bam ]; then
echo ${OUTFILE}.bam exists: skipping... 
else
echo Saving alignment to $OUTFILE ..
bwa-mem2 mem -t $NTHREAD ${BASEDIR}/02_reference/Ref_plus_biselli.fasta $READ1 $READ2 | samtools view -b -F 4 -o ${OUTFILE}.bam -
fi
done

#Count all alignments (Qual >0) on Flavobacterium fasta
for FILE in ${ALDIR}/*.bam
do
OUTFILE=${FILE/.bam/_highqual.txt}
samtools view $FILE | awk '$5>0' | cut -f3 | sort | uniq -c > $OUTFILE
OUTFILE=${FILE/.bam/_all.txt}
samtools view $FILE | cut -f3 | sort | uniq -c > $OUTFILE
done

for FILE in ${ALDIR}/*highqual.txt
do
OUTFILE=${FILE/_highqual.txt/_highqual_name.txt}
Rscript ${FUNCDIR}/03_classify_mapped_refseq.r \
-M $FILE -O $OUTFILE
done


#Repeat the exercise using only the reads classified as Flavobacterium


READDIR=${BASEDIR}/06_onlyflavoreads
ALDIR=$BASEDIR/08_bwamem_onlyflavoreads
NTHREAD=2
for READ1 in ${READDIR}/*_1.fastq.gz
do
#echo $READ1
READ2=${READ1/_1.fastq.gz/_2.fastq.gz}
FNAME=$(basename $READ1)
OUTFILE=${READ1/$READDIR/$ALDIR}
OUTFILE=${OUTFILE/_L001_1.fastq.gz/}
echo Saving alignment to $OUTFILE ..
bwa-mem2 mem -t $NTHREAD ${BASEDIR}/02_reference/Ref_plus_biselli.fasta $READ1 $READ2 | samtools view -b -F 4 -o ${OUTFILE}.bam -
done

#Count all alignments (Qual >0) on Flavobacterium fasta
for FILE in ${ALDIR}/*.bam
do
OUTFILE=${FILE/.bam/_highqual.txt}
samtools view $FILE | awk '$5>0' | cut -f3 | sort | uniq -c > $OUTFILE
OUTFILE=${FILE/.bam/_all.txt}
samtools view $FILE | cut -f3 | sort | uniq -c > $OUTFILE
done

for FILE in ${ALDIR}/*highqual.txt
do
OUTFILE=${FILE/_highqual.txt/_highqual_name.txt}
Rscript ${FUNCDIR}/03_classify_mapped_refseq.r \
-M $FILE -O $OUTFILE
done

