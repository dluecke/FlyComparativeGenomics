# ReciprocalRagtag.smk workflow for reciprocal ragtag using FlyComparativeGenomics scripts

configfile: "config.yaml"

rule all:
    input:
        "RT_final/Fp_RTfinal.fa",
        "RT_final/Fh1_RTfinal.fa",
        "RT_final/Fh2_RTfinal.fa",
        "RT_final/Mp_RTfinal.fa",
        "RT_final/Mh1_RTfinal.fa",
        "RT_final/Mh2_RTfinal.fa"

rule pass1_pri_f:
    input:
        ref=config["assemblies"]["male"]["primary"],
        qry=config["assemblies"]["female"]["primary"],
        reads=config["reads"]["female"]
    output:
        "pri/female/ragtag_output-round2/ragtag.scaffold.fasta"
    shell:
        "
        mkdir -p pri/female
        cd pri/female
        ln -s {input.ref} Mpri.fa
        ln -s {input.qry} Fpri.fa
        sbatch {config[repo_location]}/Reciprocal_RagTag.slurm Mpri.fa Fpri.fa {input.reads}
        cd ../..
        "

rule pass1_pri_m:
    input:
        ref=config["assemblies"]["female"]["primary"],
        qry=config["assemblies"]["male"]["primary"],
        reads=config["reads"]["male"]
    output:
        "pri/male/ragtag_output-round2/ragtag.scaffold.fasta"
    shell:
        "
        mkdir -p pri/male
        cd pri/male
        ln -s {input.ref} Fpri.fa
        ln -s {input.qry} Mpri.fa
        sbatch {config[repo_location]}/Reciprocal_RagTag.slurm Fpri.fa Mpri.fa {input.reads}
        cd ../..
        "

rule pass1_hapsF_hap1:
    input:
        ref=config["assemblies"]["female"]["hap2"],
        qry=config["assemblies"]["female"]["hap1"],
        reads=config["reads"]["female"]
    output:
        "haps_female/hap1/ragtag_output-round2/ragtag.scaffold.fasta"
    shell:
        "
        mkdir -p haps_female/hap1
        cd haps_female/hap1
        ln -s {input.ref} Fh2.fa
        ln -s {input.qry} Fh1.fa
        sbatch {config[repo_location]}/Reciprocal_RagTag.slurm Fh1.fa Fh2.fa {input.reads}
        cd ../..
        "

rule pass1_hapsF_hap2:
    input:
        ref=config["assemblies"]["female"]["hap1"],
        qry=config["assemblies"]["female"]["hap2"],
        reads=config["reads"]["female"]
    output:
        "haps_female/hap2/ragtag_output-round2/ragtag.scaffold.fasta"
    shell:
        "
        mkdir -p haps_female/hap2
        cd haps_female/hap2
        ln -s {input.ref} Fh1.fa
        ln -s {input.qry} Fh2.fa
        sbatch {config[repo_location]}/Reciprocal_RagTag.slurm Fh2.fa Fh1.fa {input.reads}
        cd ../..
        "

rule pass1_hapsM_hap1:
    input:
        ref=config["assemblies"]["male"]["hap2"],
        qry=config["assemblies"]["male"]["hap1"],
        reads=config["reads"]["male"]
    output:
        "haps_male/hap1/ragtag_output-round2/ragtag.scaffold.fasta"
    shell:
        "
        mkdir -p haps_male/hap1
        cd haps_male/hap1
        ln -s {input.ref} Fh2.fa
        ln -s {input.qry} Fh1.fa
        sbatch {config[repo_location]}/Reciprocal_RagTag.slurm Fh1.fa Fh2.fa {input.reads}
        cd ../..
        "

rule pass1_hapsM_hap2:
    input:
        ref=config["assemblies"]["male"]["hap1"],
        qry=config["assemblies"]["male"]["hap2"],
        reads=config["reads"]["male"]
    output:
        "haps_male/hap2/ragtag_output-round2/ragtag.scaffold.fasta"
    shell:
        "
        mkdir -p haps_male/hap2
        cd haps_male/hap2
        ln -s {input.ref} Fh1.fa
        ln -s {input.qry} Fh2.fa
        sbatch {config[repo_location]}/Reciprocal_RagTag.slurm Fh2.fa Fh1.fa {input.reads}
        cd ../..
        "

rule fake_end:
    input:
        "pri/female/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri/male/ragtag_output-round2/ragtag.scaffold.fasta",
        "haps_female/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        "haps_female/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        "haps_male/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        "haps_male/hap2/ragtag_output-round2/ragtag.scaffold.fasta"
    output:
        "RT_final/Fp_RTfinal.fa",
        "RT_final/Fh1_RTfinal.fa",
        "RT_final/Fh2_RTfinal.fa",
        "RT_final/Mp_RTfinal.fa",
        "RT_final/Mh1_RTfinal.fa",
        "RT_final/Mh2_RTfinal.fa"
    shell:
        "mkdir RT_final"