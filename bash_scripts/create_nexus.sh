#!/bin/bash

output_file="complete.nxs"

rm $output_file #just cleaning out previous iteration if running multiple times
cp -r mb mb.backup
rm -r mb
touch $output_file
echo "#NEXUS" >> $output_file
echo "BEGIN DATA;" >> $output_file
sed -n 1p dimensions >> $output_file
echo "  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=YES;
  MATRIX" >> $output_file 

# Loop through each .nxs file in the nexus_files directory
for nexus_file in nexus_files/*.nxs; do
    # Extract the gene name from the file name
    gene_name=$(basename "$nexus_file" .nxs)
    
    # Add the header for the gene
    echo "   [$gene_name]" >> $output_file
    
    # Extract the lines starting from line 7 and excluding the last two lines
    total_lines=$(wc -l < "$nexus_file")
    end_line=$((total_lines - 3))
    sed -n "7,${end_line}p" "$nexus_file" >> $output_file
    
    # Add a newline for separation (optional)
    echo "" >> $output_file
done

echo "  ;
END;

begin mrbayes;" >> $output_file
echo " " >> $output_file
sed -n 1,13p charset >> $output_file
echo " " >> $output_file
sed -n 1p partitions >> $output_file
echo "	    set partition = by_gene;
" >> $output_file
echo "	    outgroup Vampyroteuthis_infernalis;
	    " >> $output_file
sed -n 1p nucmodel >> $output_file
echo "	    prset ratepr=variable;
	    unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);" >> $output_file

cat models_to_write >> $output_file

echo "
      mcmcp ngen=100000000 printfreq=1000 samplefreq=100 stoprule=yes stopval=0.001 nruns=5 nchains=2 burninfrac=0.25;
      mcmc;
      sump;
      sumt filename=mb/$output_file contype=allcompat conformat=simple;

end;
" >> $output_file

mkdir mb
cp $output_file mb/$output_file
echo "Nexus file created."