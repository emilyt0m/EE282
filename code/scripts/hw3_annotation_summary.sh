srun -A class-ee282 --pty bash -i
conda activate ee282
wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/dmel-all-r6.48.gtf.gz
wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/md5sum.txt

md5sum -c md5sum.txt

bioawk -c gff '{print $feature}' dmel-all-r6.48.gtf.gz | sort -k3,3 | uniq -c | sort -nrk1,1

bioawk -c gff '{print $seqname, $feature}' dmel-all-r6.48.gtf.gz | grep '\<gene' | grep '[234LRXY]' | grep -v '[\d{15}]' | sort | uniq -c | sort -nrk1,1 