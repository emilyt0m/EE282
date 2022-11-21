srun -A class-ee282 --pty bash -i
conda activate ee282
wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz
wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/md5sum.txt
md5sum -c <(grep all-chrom md5sum.txt)

faSize dmel-all-chromosome-r6.48.fasta.gz > dmel-all-chromosome-summary.txt

cat dmel-all-chromosome-summary.txt
