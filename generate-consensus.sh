me=`basename "$0"`

command -v bwa >/dev/null 2>&1 || { echo "I require bwa but it's not installed.  Aborting." >&2; }
command -v samtools >/dev/null 2>&1 || { echo "I require samtools but it's not installed.  Aborting." >&2; }
command -v bcftools >/dev/null 2>&1 || { echo "I require bcftools but it's not installed.  Aborting." >&2; }
command -v vcfutils.pl >/dev/null 2>&1 || { echo "I require vcfutils.pl but it's not installed.  Aborting." >&2; }

if [ $# -ne 4 ]; then
    echo Script needs directory input.
    echo Script usage: $me ./path/to/contigs ./path/to/output-dir ./path/to/reference_fasta_file.fasta num_processors
    exit 1
fi

workdir=$(readlink -e  $1)
outputdir=$(readlink -f $2)
reference=$(readlink -e $3)
proc_num=$4

mkdir $outputdir
 
REGEX="*.contigs.fasta" 

bwa index $reference

for ARQ in $workdir/*
do
    if [[ $ARQ == $REGEX ]]; then
        completename=`basename "$ARQ"`;
        IFS='.' read -r -a array <<< "$completename"
        samp=${array[0]}
        echo "Working on :" $completename

        bwa mem -t $proc_num $reference $ARQ | samtools sort -@$proc_num -O BAM -o $outputdir/$samp.bam - 
        samtools mpileup -E -uf $reference $outputdir/$samp.bam | bcftools call -c | vcfutils.pl vcf2fq > $outputdir/$samp.fastq
        
    fi	
done




