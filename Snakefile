configfile: "config.yaml"
workdir: "/scif/data"

SAMPLES, = glob_wildcards("raw_data/{sample}_R1.fastq.gz")

rule all:
    input:
        expand("tag_directories/{sample}/track_info.txt",sample=SAMPLES),
        expand("tag_directories/{sample}/peaks.txt",sample=SAMPLES)

rule bowtie2_map:
    input:
        "raw_data/{sample}_R1.fastq.gz"
    output:
        "aligned_files/{sample}.sam"
    params:
        "bowtie_index/idx"
    threads:
        50
    shell:
        "scif run bowtie2 '-p {threads} -x $SCIF_DATA/{params} $SCIF_DATA/{input} > $SCIF_DATA/{output}'"

rule sort_sam:
    input:
        "aligned_files/{sample}.sam"
    output:
        "aligned_files/{sample}.sorted.bam"
    threads:
        50
    shell:
        """
        scif run samtools 'sort -o $SCIF_DATA/{output} -@ {threads} $SCIF_DATA/{input}'
        rm {input}
        """

rule make_tag_directory:
    input:
        "aligned_files/{sample}.sorted.bam"
    output:
        "tag_directories/{sample}"
    shell:
        "scif run makeTagDirectory '$SCIF_DATA/{output} $SCIF_DATA/{input}'"

rule make_bigwig:
    input:
        "tag_directories/{sample}"
    output:
        "tag_directories/{sample}/track_info.txt"
    params:
        "chrom.sizes"
    shell:
        "scif run makeUCSCfile '$SCIF_DATA/{input} -o auto -bigWig $SCIF_DATA/{params} -fsize 1e20 > $SCIF_DATA/{output}'"

rule find_peaks:
    input:
        "tag_directories/{sample}"
    output:
        "tag_directories/{sample}/peaks.txt"
    params:
        config["peak_cmds"]
    shell:
        "scif run findPeaks '$SCIF_DATA/{input} -o $SCIF_DATA/{output} {params}'"