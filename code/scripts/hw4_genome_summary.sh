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

library(ggplot2)
length_short <- read.table("length_short.txt", header = FALSE)
length_long <- read.table("length_long.txt", header = FALSE)
colnames(length_short) <- c('name', 'length', 'assembly' )
colnames(length_long) <- c('name', 'length', 'assembly' )

plot_length_short <- ggplot(data = length_short, aes (x=log_length)) + geom_histogram(bins=200) + labs(title = 'Log Sequence Length, <100kb')plot_length_short

plot_length_long <- ggplot(data = length_long, aes (x=log_length)) + geom_histogram(bins=7) + labs(title = 'Log Sequence Length, >=100kb')plot_length_long

bioawk -c fastx ' { print length($seq) } ' less.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > seq_less.lengths
plotCDF2 seq_less.lengths CDFless.png

bioawk -c fastx ' { print length($seq) } ' more.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > seq_more.lengths
plotCDF2 seq_more.lengths CDFmore.png

minimap2 -x ava-ont -t16 /pub/etom2/classrepos/EE282/data/raw/iso1_onp_a2_1kb.fastq{,} | gzip -1 > /pub/etom2/classrepos/EE282/data/raw/iso1_onp_reads.paf.gz &

miniasm -f /pub/etom2/classrepos/EE282/data/raw/iso1_onp_a2_1kb.fastq /pub/etom2/classrepos/EE282/data/raw/iso1_onp_reads.paf.gz > /pub/etom2/classrepos/EE282/data/processed/iso1_assembly.gfa

n50 () {
  bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' $1 \
  | sort -rn \
  | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
}
awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' iso1_assembly.gfa \
| tee >(n50 /dev/stdin > n50.txt) \
| fold -w 60 \
> unitigs.fa

r6url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/215/GCA_000001215.4_Release_6_plus_ISO1_MT/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz"

wget -O - -q $r6url \
| tee >( \
  bioawk -c fastx ' { print length($seq) } ' \
  | sort -rn \
  | awk ' BEGIN { print "Assembly\tLength\nFB_Scaff\t0" } { print "FB_Scaff\t" $1 } ' \
  > data/processed/ISO1.r6.scaff.sorted.sizes.txt & ) \
| faSplitByN /dev/stdin /dev/stdout 10 \
| bioawk -c fastx ' { print length($seq) } ' \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nFB_Ctg\t0" } { print "FB_Ctg\t" $1 } ' \
> data/processed/ISO1.r6.ctg.sorted.sizes.txt &

bioawk -c fastx ' { print length($seq) } ' unitigs.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nMinimap_Ctg\t0" } { print "Minimap_Ctg\t" $1 } ' \
> data/processed/minimap.ctg.sorted.sizes.txt

plotCDF2 data/processed/*.sizes.txt output/figures/CDFcomparison.png

busco -c 16 -i unitigs.fa -l diptera_odb10 -o Dmel_busco -m genome

busco -c 32 -i ISO1.r6.ctg.fa -l diptera_odb10 -o DmelNCBI_busco -m genome
