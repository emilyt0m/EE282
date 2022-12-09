srun -A class-ee282 --pty bash -i
conda activate ee282
wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz
gunzip dmel-all-chromosome-r6.48.fasta
faFilter -maxSize=100000 dmel-all-chromosome-r6.48.fasta less.fasta
faSize less.fasta
faFilter -minSize=100001 dmel-all-chromosome-r6.48.fasta more.fasta
faSize more.fasta


bioawk -c fastx '{ print $name "\t" length($seq) "\tshort" } ' less.fasta > length_short.txt
bioawk -c fastx '{ print $name "\t" length($seq) "\tlong" } ' more.fasta > length_long.txt

bioawk -c fastx '{ print $name "\t" gc($seq) "\tshort" } ' less.fasta > gc_short.txt
bioawk -c fastx '{ print $name "\t" gc($seq) "\tlong" } ' more.fasta > gc_long.txt

bioawk -c fastx ' { print length($seq) } ' less.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > seq_less.lengths
plotCDF2 seq_less.lengths CDFless.png

bioawk -c fastx ' { print length($seq) } ' more.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > seq_more.lengths
plotCDF2 seq_more.lengths CDFmore.png

