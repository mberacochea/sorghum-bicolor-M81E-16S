#!/bin/bash

# params
SEQ=$1
MANIFEST=$2

# helpers
run()
{
  echo -e \\n"$@"\\n
  eval "$@"
}

print()
{
    echo -e \\n"${REV}${1}${NORM}"\\n
}

create_folder()
{
    if [ ! -d $1 ]; then
		mkdir $1
    fi
}

delete_prev_results()
{
    if [ -f $1 ]; then
		rm $1
    fi
}

usage()
{
	echo -e "Usage: $(basename "$0") <samples.fastq> <manifest.txt>\n"
	exit 1
}

if [ -z $1 ];then
	print "Provide a fastq file."
	usage
fi

if [ -z $2 ];then
	print "Provide a Manifest file."
	usage
fi

# fonts
NORM=`tput sgr0`
REV=`tput smso`

# output
FOLDER="res_"$SEQ"_$(date +'%d-%m-%Y')"

create_folder $FOLDER

RES_PWD=$PWD/$FOLDER

if [ ! -d $RES_PWD/fastas ]; then
    # create dir
    mkdir $RES_PWD/fastas
    # copy fastq
    cp $PWD/$SEQ.fasta $RES_PWD/fastas/$SEQ.fasta
fi

print "Working with ${1}"
print "BMP workflow http://www.brmicrobiome.org/#!16s-profiling-ion-torrent-new/cpdg"

#-----------------------------------------------------------------------------------------------------#
print "Dereplication <<<USING USEARCH 8>>>"

R_DEREP="${SEQ}_derep.fa"

# usearch
run "usearch8 -derep_fulllength ${RES_PWD}/fastas/${SEQ}.fasta -fastaout ${RES_PWD}/fastas/${R_DEREP} -sizeout"

#-----------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------#

print "Abundance sort and discard singletons <<<USING USEARCH 8>>>"

R_ABUND="${SEQ}_sorted.fa"
R_ABUND_FULL=$RES_PWD/fastas/$R_ABUND

delete_prev_results $R_ABUND_FULL

run "usearch8 -sortbysize $RES_PWD/fastas/$OUT_3 -fastaout $R_ABUND_FULL -minsize 2"

#-----------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------#

print "OTU clustering using UPARSE method <<<USING USEARCH 8>>>"

R_OTU_CLUS="${SEQ}_otus1.fa"
R_OTU_CLUS_DIR="${RES_PWD}/bmp_otus"
R_OTU_CLUS_FULL=$R_OTU_CLUS_DIR/$R_OTU_CLUS

# delete previuos runs content
create_folder $R_OTU_CLUS_DIR

delete_prev_results $R_OTU_CLUS_FULL

# usearch
run "usearch8 -cluster_otus $R_ABUND_FULL -otus $R_OTU_CLUS_FULL"

#-----------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------#

print "Chimera filtering using reference database <<<USING USEARCH 8>>>"

R_CHIMERA_F="${SEQ}_otus.fa"
R_CHIMERA_F_DIR=$R_OTU_CLUS_DIR
R_CHIMERA_F_FULL=$R_CHIMERA_F_DIR/$R_CHIMERA_F

# copy the rdp_gold.fa file
if [ -f $PWD/rdp_gold.fa ]; then
    echo -e "\nCopying rdp_gold.fa"
    cp $PWD/rdp_gold.fa $RES_PWD/rdp_gold.fa
fi

# copy
if [ -f ../$PWD/rdp_gold.fa ]; then
    echo -e "\nCopying rdp_gold.fa"
    cp ../$PWD/rdp_gold.fa $RES_PWD/rdp_gold.fa
fi

if [ ! -f $RES_PWD/rdp_gold.fa ]; then
    echo -e "\nDownloading rdg_gold.fa"
    wget "http://drive5.com/uchime/rdp_gold.fa" -O $RES_PWD/rdp_gold.fa
fi

delete_prev_results $R_CHIMERA_F_FULL

# usearch
run "usearch8 -uchime_ref ${R_OTU_CLUS_FULL} -db ${RES_PWD}/rdp_gold.fa -strand plus -minh 1.0 -nonchimeras ${R_CHIMERA_F_FULL}"

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#

print "Map reads back to OTU database <<<USING USEARCH8>>>"

R_OTU_MAP="${SEQ}_map.biom"
R_OTU_MAP_DIR="${RES_PWD}/bmp_map"
R_OTU_MAP_FULL=$R_OTU_MAP_DIR/$R_OTU_MAP

create_folder $R_OTU_MAP_DIR

delete_prev_results $R_OTU_MAP_DIR/$R_OTU_MAP

# usearch
run "usearch8 -usearch_global ${RES_PWD}/fastas/${SEQ}.fasta -db ${R_CHIMERA_F_FULL} -strand plus -id 0.97 -biomout ${R_OTU_MAP_FULL}"

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#

print "Assign taxonomy to OTUS using uclust method on QIIME (use the file “otus.fa” from UPARSE as input file)"

R_TAX="${SEQ}_tax"
R_TAX_DIR=$RES_PWD/"tax"
R_TAX_FULL=$R_TAX_DIR/$R_TAX

create_folder $R_TAX_DIR

delete_prev_results $R_TAX_FULL

# qiime
run "assign_taxonomy.py -i ${R_CHIMERA_F_FULL} -o ${R_TAX_FULL}"

#----------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------#

print "Align sequences on QIIME, using greengenes reference sequences (use the file “otus.fa” from UPARSE as input file)"

R_ALIGN="${SEQ}_align"
R_ALIGN_DIR=$RES_PWD/"aligned"
R_ALIGN_FULL=$R_ALIGN_DIR/$R_ALIGN

create_folder $R_ALIGN_DIR

delete_prev_results $R_ALIGN_FULL

# qiime
run "align_seqs.py -i ${R_CHIMERA_F_FULL} -o ${R_ALIGN_FULL}"

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#

print "Filter alignments on QIIME"

R_ALIGN_F="${SEQ}_filt_align"
R_ALIGN_F_DIR=$RES_PWD/"filt_aligned"
R_ALIGN_F_FULL=$R_ALIGN_F_DIR/$R_ALIGN_F

create_folder $R_ALIGN_F_DIR

delete_prev_results $R_ALIGN_F_FULL

# qiime
run "filter_alignment.py -i ${R_ALIGN_FULL}/${SEQ}_otus_aligned.fasta -o ${R_ALIGN_F_FULL}"

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#

print "Make the reference tree on QIIME"

R_TREE="${SEQ}_tre"
R_TREE_DIR=$RES_PWD/"phylogeny"
R_TREE_FULL=$R_TREE_DIR/$R_TREE

create_folder $R_TREE_DIR

delete_prev_results $R_TREE_FULL


# qiime
run "make_phylogeny.py -i ${R_ALIGN_F_FULL}/${SEQ}_otus_aligned_pfiltered.fasta -o ${R_TREE_FULL}"

print "Add metadata (taxonomy) to OTU table"

R_TAX_METADATA="${SEQ}_otu_table_tax.biom"

# biom
run "biom add-metadata -i $R_OTU_MAP_FULL -o ${RES_PWD}/${R_TAX_METADATA} --observation-metadata-fp ${R_TAX_FULL}/${SEQ}_otus_tax_assignments.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy --float-fields confidence"

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#

print "Check OTU Table on QIIME."

R_OTU_TABLE="${SEQ}_results_biom_table_summary.txt"

# biom
run "biom summarize-table -i ${RES_PWD}/${R_TAX_METADATA} -o ${RES_PWD}/${R_OTU_TABLE}"

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------#

print " Run diversity analyses on QIIME (or any other analysis of your choice). The parameter '-e' is the sequencing depth to use for even sub-sampling and maximum rarefaction depth. You should review the output of the ‘biom summarize-table’ command to decide on this value."

# ask the user for sampling depth
VALID=0

echo -e "\nSelected manifest: "$MANIFEST

echo -e "\nSampling depth"

VALID=0
E=0
while [ $VALID == 0 ]; do
    read e
    if [ ! -z "${e}" ]; then
	VALID=1
    else
	echo "e is not valid" $PWD
    fi
done
E=$e

# parse the summary and get the -e value

# qiime
run "core_diversity_analyses.py -i ${RES_PWD}/${R_TAX_METADATA} -m ${PWD}/${MANIFEST} -t ${R_TREE_FULL} -o ${RES_PWD}/diversity_analysis -e ${E} --suppress_beta_diversity"
