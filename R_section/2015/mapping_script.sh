for file in mapped/*
do
	if [ ${file: -4} == ".sam" ]; then # check that it does indeed end in the .sam file extension
		filename=$(basename $file) # omit the directory leading up the file name
		filenamenoext="${filename%.*}" # Get the filename without an extension
		# create a new file with the same name but in directory mapped_fixed/
		# directory was already created through "mkdir mapped_fixed"
	
		echo "Working on" $filename"..." 
		awk 'BEGIN {OFS="\t"}; {gsub("^.*NC_", "|NC_", $3); gsub("\\|", "", $3); print}' mapped/$filename > mapped_fixed/$filename
		# It isn't as good as I wanted... It would have been nice to not include "|NC_" as part of the subtitution
		# Instead of replace with itself. Would have been better to be able to select everything before that.
		# But testing this on ERR459206.sam appeared to do the job...

		# Now the counts. Do htseq-count on the mapped_fixed
		python -m HTSeq.scripts.count -f sam -t gene -i locus_tag -a 0 mapped_fixed/$filename mapped/yeast.gff > htseqcounts/$filenamenoext.count.txt

		# Cleanup. Delete the mapped_fixed sam file
		rm mapped_fixed/$filename

		echo "Done" $filename"."
	fi
done

# Cleanup. Delete the now-empty mapped_fixed directory
rmdir mapped_fixed
