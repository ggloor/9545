BIN="/Groups/twntyfr/bin/bowtie2-2.1.0/"

# Build index
# $BIN/bowtie2-build complete_genome.fna index/complete

for file in fastq_files/*
do
	# Get the name of the file only without the folders, but still has extension
	filename=$(basename $file)
	# Get rid of the .gz
	filename="${filename%.*}"
	# Get rid of the .fastq
	filename="${filename%.*}"
	# Execute the command that's the whole point of this
	# Edited "$BIN/bowtie2" to "$BIN"bowtie2"" because the original was writing ...//bowtie2. Not sure if that matters
	echo $filename	
	$BIN/bowtie2 -x yeast/index/complete -U fastq_files/$filename.fastq.gz -p 8 -S mapped/$filename.sam
done
