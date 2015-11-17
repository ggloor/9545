ID="ERR458"
NUM=990
FA=".fastq.gz"

while [ $NUM -lt 1000 ]; do 
	echo "$ID$NUM/$ID$NUM$FA"
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$ID/$ID$NUM/$ID$NUM$FA
	let NUM=NUM+1
done
