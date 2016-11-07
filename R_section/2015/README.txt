# bin folder is:
BIN="/Groups/twntyfr/bin/bowtie2-2.1.0/"

# mkdir index
# then build the index with bowtie2
$BIN/bowtie2-build complete_genome.fna index/complete

# example mapping command:
barton_RNAseq ggloor$ $BIN/bowtie2 -x yeast/index/complete -U fastq_files/ERR459000.fastq.gz -p 8 -S mapped/ERR459000.sam

exit 1

# FOR REFERENCE ONLY - this is Jean's shell script for mapping to corn
#Map all files to maize genome

#nohup /Volumes/bin/bowtie2-2.2.6/bowtie2 -x ../reference_sequences/GCF_000005005 -U ~/Downloads/GT1_A_R1.fastq.gz -p 18 -S GT1_o
ut.sam &

BIN=/Volumes/bin
DATA_DIR=/home/mmacklai/Downloads
WORKING_DIR=/Volumes/data/A_n_L/sept2015

echo "# Mapper version:"
$BIN/bowtie2 --version
echo -e "\n"

for f in $DATA_DIR/*_R1.fastq.gz; do
#	echo "$f"
	B=`basename $f`
#	echo "basename: $B"
	NAME=`echo $B | cut -d "_" -f1`
#	echo "name: $NAME"

	echo "# Running $NAME"
	$BIN/bowtie2 -x /Volumes/data/A_n_L/reference_sequences/GCF_000005005 -U $f -p 30 -S $WORKING_DIR/${NAME}_maize.sam
	echo -e "\n"

done

