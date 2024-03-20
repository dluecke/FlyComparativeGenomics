rule extract_hifi_reads:
    input:
        "raw_hifi/{sample}.bam"
    output:
        fa="hifi_fasta/{sample}.fasta",
        fq="hifi_fastq/{sample}.fastq"
    shell:
        "module load samtools/1.17; "
        "samtools fasta {input} > {output.fa}; "
        "samtools fastq {input} > {output.fq}"

rule fcs_adaptor:
    input:
        "hifi_fasta/{sample}.fasta"
    output:
        "fcs-adaptor/{sample}/fcs_adaptor_report.txt"
    shell:
        "/project/vpgru/software/fcsadaptor/run_fcsadaptor.sh "
        "--fasta-input {input} "
        "--output-dir ./fcs-adaptor/{sample}/ "
        "--euk --container-engine singularity "
        "--image /project/vpgru/software/fcsadaptor/fcs-adaptor.sif"

rule HiFiAdapterFilt:
    input:
        fq="hifi_fastq/{sample}.fastq", 
        fcs="fcs-adaptor/{sample}/fcs_adaptor_report.txt"
    output:
        "filtered_hifi/{sample}.filt.fastq.gz"
    shell:
        "module load bamtools/2.5.2; "
        "module load pigz/2.7; "
        "/project/vpgru/software/ HiFiAdapterFilt/hifiadapterfiltFCS.sh "
        "-f {output.fcs} -r {output.fq} -o filtered_hifi"



