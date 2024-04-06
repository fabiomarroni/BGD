#!/bin/bash

#Variables you need to change to your local setting
###################################################
FUNCDIR=/full/path/to/functions/
BASEDIR=/base/directory/
###################################################



############################################
#
#  _              _              ___  
# | |            | |            |__ \ 
# | | ___ __ __ _| | _____ _ __    ) |
# | |/ / '__/ _` | |/ / _ \ '_ \  / / 
# |   <| | | (_| |   <  __/ | | |/ /_ 
# |_|\_\_|  \__,_|_|\_\___|_| |_|____|
#                                     
############################################
                                     

#We assume you already saved the fastq files in READDIR folder
READDIR=${BASEDIR}/00_reads
KRES=${BASEDIR}/05_kraken2

######################
#
# kraken2 on refseq
#
######################



#Create a kraken database from the 16s refseq database
KRES=${BASEDIR}/06_kraken2_16s_RefSeq
KDB=/projects/marroni/databases/16S_RefSeq_k2db
mkdir -p $KRES

#Download taxonomy
kraken2-build --download-taxonomy --db $KDB

#Add fasta to library
kraken2-build --add-to-library $KDB/RefSeq_16s.fasta --db $KDB

#Actually build the database
kraken2-build --build --db $KDB

READPATTERN=_R1_001.fastq.gz
#End of parameters that you likely need to change

#Run kraken 
for aaa in $READDIR/*${READPATTERN}
do
read1=$(basename $aaa)
read2=${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}
pref=${read1/_R1_001.fastq.gz/}
echo $read1
echo $read2
mkdir -p ${KRES}/logs
if [ -s ${KRES}/${pref}.kraken.report.txt ]
	then
	echo "Kraken already run, i will skip it"	
else
cd ${KRES} 
export TMPDIR=${KRES}
kraken2 --threads 12 --paired --gzip-compressed --db $KDB $READDIR/${read1} $READDIR/${read2} --output ${KRES}/${pref}.kraken --classified-out ${KRES}/${pref}#.fastq --use-names --report ${KRES}/${pref}.kraken.report.txt
fi
done


###########################################                                      
#  _                    _              
# | |                  | |             
# | |__  _ __ __ _  ___| | _____ _ __  
# | '_ \| '__/ _` |/ __| |/ / _ \ '_ \ 
# | |_) | | | (_| | (__|   <  __/ | | |
# |_.__/|_|  \__,_|\___|_|\_\___|_| |_|
#                                      
###########################################                                      

#Build bracken DB for 16s RefSeq
bracken-build -d $KDB -t 12 -k 35 -l 250 


THRESHOLD=10
READ_LEN=250
KRES=${BASEDIR}/06_kraken2_16s_RefSeq
for CLASSIFICATION_LEVEL in S G
do
for READ1 in ${READDIR}/*${READPATTERN}
do
read1=$(basename $READ1)
read1=${read1/_R1_001.fastq.gz/}
echo $read1
#Run bracken on kraken results
echo $THRESHOLD
bracken -d ${KDB} -i ${KRES}/${read1}.kraken.report.txt -o ${KRES}/${read1}_${CLASSIFICATION_LEVEL}.bracken.txt -r ${READ_LEN} -l ${CLASSIFICATION_LEVEL} -t ${THRESHOLD}
done
done




################################################
#  _____  _       _           _                     _                      
# |  __ \| |     | |         | |                   | |                     
# | |__) | | ___ | |_    __ _| |__  _   _ _ __   __| | __ _ _ __   ___ ___ 
# |  ___/| |/ _ \| __|  / _` | '_ \| | | | '_ \ / _` |/ _` | '_ \ / __/ _ \
# | |    | | (_) | |_  | (_| | |_) | |_| | | | | (_| | (_| | | | | (_|  __/
# |_|    |_|\___/ \__|  \__,_|_.__/ \__,_|_| |_|\__,_|\__,_|_| |_|\___\___|
#                                                                          
################################################
                                                                          
#Using minikraken db
Rscript ${FUNCDIR}/01_summarize_bracken.r \
       -N 15 --removeme '' \
--outgraph ${BASEDIR}/03_plots/species_bracken_barplot.pdf

Rscript ${FUNCDIR}/01_summarize_bracken.r \
       -N 15 --removeme 10_SANO \
	   --outgraph ${BASEDIR}/03_plots/species_bracken_clean_barplot.pdf

#Using RefSeq16s db
Rscript ${FUNCDIR}/01_summarize_bracken.r \
    --indir ${BASEDIR}/06_kraken2_16s_RefSeq/ \
    -N 15 --removeme '' \
    --outfile ${BASEDIR}/04_tables_16s_RefSeq/bracken.txt \
    --outgraph ${BASEDIR}/03_plots_16s_RefSeq/species_bracken_barplot.pdf

Rscript ${FUNCDIR}/01_summarize_bracken.r \
    --indir ${BASEDIR}/06_kraken2_16s_RefSeq/ \
    -N 15 --removeme 10_SANO \
    --outfile ${BASEDIR}/04_tables_16s_RefSeq/bracken.txt \
    --outgraph ${BASEDIR}/03_plots_16s_RefSeq/species_bracken_clean_barplot.pdf


#############
#
# Genus level
#
#############
#Using minikraken db
Rscript ${FUNCDIR}/01a_summarize_genus_bracken.r \
    --indir /projects/marroni/BGD/05_kraken2 \
    --outfile /projects/marroni/BGD/04_tables/genus_bracken.txt \
    -N 15 --removeme '' \
    --outgraph ${BASEDIR}/03_plots/genus_bracken_barplot.pdf

Rscript ${FUNCDIR}/01a_summarize_genus_bracken.r \
    --indir /projects/marroni/BGD/05_kraken2 \
    --outfile /projects/marroni/BGD/04_tables/genus_bracken.txt \
       -N 15 --removeme 10_SANO \
	   --outgraph ${BASEDIR}/03_plots/genus_bracken_clean_barplot.pdf

#Using RefSeq16s db
Rscript ${FUNCDIR}/01a_summarize_genus_bracken.r \
    --indir ${BASEDIR}/06_kraken2_16s_RefSeq/ \
    -N 15 --removeme '' \
    --outfile ${BASEDIR}/04_tables_16s_RefSeq/genus_bracken.txt \
    --outgraph ${BASEDIR}/03_plots_16s_RefSeq/genus_bracken_barplot.pdf

Rscript ${FUNCDIR}/01a_summarize_genus_bracken.r \
    --indir ${BASEDIR}/06_kraken2_16s_RefSeq/ \
    -N 15 --removeme 10_SANO \
    --outfile ${BASEDIR}/04_tables_16s_RefSeq/genus_bracken.txt \
    --outgraph ${BASEDIR}/03_plots_16s_RefSeq/genus_bracken_clean_barplot.pdf



###############################################################
#  _____               _____  ______  _____            ___  
# |  __ \             |  __ \|  ____|/ ____|          |__ \ 
# | |__) |   _ _ __   | |  | | |__  | (___   ___  __ _   ) |
# |  _  / | | | '_ \  | |  | |  __|  \___ \ / _ \/ _` | / / 
# | | \ \ |_| | | | | | |__| | |____ ____) |  __/ (_| |/ /_ 
# |_|  \_\__,_|_| |_| |_____/|______|_____/ \___|\__, |____|
#                                                   | |     
#                                                   |_|     
#
##################################################################

############################
#
# Genus level
#
############################


#All samples on 16s_RefSeq
for COND in Status Stage StageStatus
do
mkdir -p ${BASEDIR}/07_DA_16s_RefSeq_genus_${COND}
Rscript ${FUNCDIR}/02_DEseq.r \
--abundance ${BASEDIR}/04_tables_16s_RefSeq/genus_bracken_raw.txt \
--removeme '' \
-C ${BASEDIR}/01_data/ID2400_16S_4th_analysis_config_u.tsv \
-c $COND -r 200 -O ${BASEDIR}/07_DA_16s_RefSeq_genus_${COND}/
done

#Only pre
for COND in Status 
do
mkdir -p ${BASEDIR}/07_DA_16s_RefSeq_onlypre_genus_${COND}
Rscript ${FUNCDIR}/02_DEseq.r \
--abundance ${BASEDIR}/04_tables_16s_RefSeq/genus_bracken_raw.txt \
--removeme '' \
-C ${BASEDIR}/01_data/ID2400_16S_4th_analysis_config_u_only_pre.tsv \
-c $COND -r 200 -O ${BASEDIR}/07_DA_16s_RefSeq_onlypre_genus_${COND}/
done

#Only post
for COND in Status 
do
mkdir -p ${BASEDIR}/07_DA_16s_RefSeq_onlypost_genus_${COND}
Rscript ${FUNCDIR}/02_DEseq.r \
--abundance ${BASEDIR}/04_tables_16s_RefSeq/genus_bracken_raw.txt \
--removeme '' \
-C ${BASEDIR}/01_data/ID2400_16S_4th_analysis_config_u_only_post.tsv \
-c $COND -r 200 -O ${BASEDIR}/07_DA_16s_RefSeq_onlypost_genus_${COND}/
done

#No outlier, on 16s_RefSeq
for COND in Status Stage StageStatus
do
mkdir -p ${BASEDIR}/07_DA_clean_16s_RefSeq_genus_${COND}
Rscript ${FUNCDIR}/02_DEseq.r \
--abundance ${BASEDIR}/04_tables_16s_RefSeq/genus_bracken_raw.txt \
--removeme '10' \
-C ${BASEDIR}/01_data/ID2400_16S_4th_analysis_config_u.tsv \
-c $COND -r 200 -O ${BASEDIR}/07_DA_clean_16s_RefSeq_genus_${COND}/
done

for COND in Status 
do
mkdir -p ${BASEDIR}/07_DA_clean_16s_RefSeq_onlypre_genus_${COND}
Rscript ${FUNCDIR}/02_DEseq.r \
--abundance ${BASEDIR}/04_tables_16s_RefSeq/genus_bracken_raw.txt \
--removeme '10' \
-C ${BASEDIR}/01_data/ID2400_16S_4th_analysis_config_u_only_pre.tsv \
-c $COND -r 200 -O ${BASEDIR}/07_DA_clean_16s_RefSeq_onlypre_genus_${COND}/
done

for COND in Status 
do
mkdir -p ${BASEDIR}/07_DA_clean_16s_RefSeq_onlypost_genus_${COND}
Rscript ${FUNCDIR}/02_DEseq.r \
--abundance ${BASEDIR}/04_tables_16s_RefSeq/genus_bracken_raw.txt \
--removeme '10' \
-C ${BASEDIR}/01_data/ID2400_16S_4th_analysis_config_u_only_post.tsv \
-c $COND -r 200 -O ${BASEDIR}/07_DA_clean_16s_RefSeq_onlypost_genus_${COND}/
done



