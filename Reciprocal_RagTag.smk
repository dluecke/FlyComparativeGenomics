# ReciprocalRagtag.smk workflow for reciprocal ragtag using FlyComparativeGenomics scripts

configfile: "config.yaml"
localrules: all, final

rule all:
    input:
        "RT_final/Fp_RTfinal.fa",
        "RT_final/Fh1_RTfinal.fa",
        "RT_final/Fh2_RTfinal.fa",
        "RT_final/Mp_RTfinal.fa",
        "RT_final/Mh1_RTfinal.fa",
        "RT_final/Mh2_RTfinal.fa"


# PASS 1: pri_F with pri_M, hap1_F with hap2_F, hap1_M with hap2_F
rule pass1_pri:
    input:
        priF=config["assemblies"]["female"]["primary"],
        priM=config["assemblies"]["male"]["primary"],
        readsF=config["reads"]["female"],
        readsM=config["reads"]["male"]
    output:
        "pri/female/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri/male/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri/female/ragtag_round2-pct_chrs.txt",
        "pri/male/ragtag_round2-pct_chrs.txt"
    shell:
        r"""
        mkdir -p pri/female pri/male
        cd pri
        ln -s {input.priF} Fpri.fa
        ln -s {input.priM} Mpri.fa
        cd female
        {config[repo_location]}/Reciprocal_RagTag.sh ../Mpri.fa ../Fpri.fa {input.readsF}
        cd ../male
        {config[repo_location]}/Reciprocal_RagTag.sh ../Fpri.fa ../Mpri.fa {input.readsM}
        cd ../..
        """

rule pass1_hapsF:
    input:
        hap1=config["assemblies"]["female"]["hap1"],
        hap2=config["assemblies"]["female"]["hap2"],
        reads=config["reads"]["female"]
    output:
        "haps_female/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        "haps_female/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        "haps_female/hap1/ragtag_round2-pct_chrs.txt",
        "haps_female/hap2/ragtag_round2-pct_chrs.txt"
    shell:
        r"""
        mkdir -p haps_female/hap1 haps_female/hap2
        cd haps_female
        ln -s {input.hap1} Fh1.fa
        ln -s {input.hap2} Fh2.fa
        cd hap1
        {config[repo_location]}/Reciprocal_RagTag.sh ../Fh2.fa ../Fh1.fa {input.reads}
        cd ../hap2
        {config[repo_location]}/Reciprocal_RagTag.sh ../Fh1.fa ../Fh2.fa {input.reads}
        cd ../..
        """

rule pass1_hapsM:
    input:
        hap1=config["assemblies"]["male"]["hap1"],
        hap2=config["assemblies"]["male"]["hap2"],
        reads=config["reads"]["male"]
    output:
        "haps_male/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        "haps_male/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        "haps_male/hap1/ragtag_round2-pct_chrs.txt",
        "haps_male/hap2/ragtag_round2-pct_chrs.txt"
    shell:
        r"""
        mkdir -p haps_male/hap1 haps_male/hap2
        cd haps_male
        ln -s {input.hap1} Mh1.fa
        ln -s {input.hap2} Mh2.fa
        cd hap1
        {config[repo_location]}/Reciprocal_RagTag.sh ../Mh2.fa ../Mh1.fa {input.reads}
        cd ../hap2
        {config[repo_location]}/Reciprocal_RagTag.sh ../Mh1.fa ../Mh2.fa {input.reads}
        cd ../..
        """


# PASS 2: output from pass 1; pri_F with hap1_F, pri_F with hap2_F, pri_M with hap1_M, pri_M with hap2_M
rule pass2_F:
    input:
        pri="pri/female/ragtag_output-round2/ragtag.scaffold.fasta",
        hap1="haps_female/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        hap2="haps_female/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        reads=config["reads"]["female"],
        pDone="pri/female/ragtag_round2-pct_chrs.txt",
        h1Done="haps_female/hap1/ragtag_round2-pct_chrs.txt",
        h2Done="haps_female/hap2/ragtag_round2-pct_chrs.txt"
    output:
        "pri_haps-female/first_hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-female/first_hap1/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-female/first_hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-female/first_hap2/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-female/first_hap1/pri/ragtag_round2-pct_chrs.txt",
        "pri_haps-female/first_hap1/hap1/ragtag_round2-pct_chrs.txt",
        "pri_haps-female/first_hap2/pri/ragtag_round2-pct_chrs.txt",
        "pri_haps-female/first_hap2/hap2/ragtag_round2-pct_chrs.txt"
    shell:
        r"""
        mkdir -p pri_haps-female/first_hap1/pri pri_haps-female/first_hap1/hap1
        mkdir -p pri_haps-female/first_hap2/pri pri_haps-female/first_hap2/hap2
        cd pri_haps-female/
        ln -s ../{input.pri} Fp_RT1Mp.fa
        cd first_hap1
        ln -s ../../{input.hap1} Fh1_RT1Fh2.fa
        cd pri
        {config[repo_location]}/Reciprocal_RagTag.sh ../Fh1_RT1Fh2.fa ../../Fp_RT1Mp.fa {input.reads}
        cd ../hap1
        {config[repo_location]}/Reciprocal_RagTag.sh ../../Fp_RT1Mp.fa ../Fh1_RT1Fh2.fa {input.reads}
        cd ../../first_hap2
        ln -s ../../{input.hap2} Fh2_RT1Fh1.fa
        cd pri
        {config[repo_location]}/Reciprocal_RagTag.sh ../Fh2_RT1Fh1.fa ../../Fp_RT1Mp.fa {input.reads}
        cd ../hap2
        {config[repo_location]}/Reciprocal_RagTag.sh ../../Fp_RT1Mp.fa ../Fh2_RT1Fh1.fa {input.reads}
        cd ../../..
        """

rule pass2_M:
    input:
        pri="pri/male/ragtag_output-round2/ragtag.scaffold.fasta",
        hap1="haps_male/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        hap2="haps_male/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        reads=config["reads"]["male"],
        pDone="pri/male/ragtag_round2-pct_chrs.txt",
        h1Done="haps_male/hap1/ragtag_round2-pct_chrs.txt",
        h2Done="haps_male/hap2/ragtag_round2-pct_chrs.txt"
    output:
        "pri_haps-male/first_hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-male/first_hap1/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-male/first_hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-male/first_hap2/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-male/first_hap1/pri/ragtag_round2-pct_chrs.txt",
        "pri_haps-male/first_hap1/hap1/ragtag_round2-pct_chrs.txt",
        "pri_haps-male/first_hap2/pri/ragtag_round2-pct_chrs.txt",
        "pri_haps-male/first_hap2/hap2/ragtag_round2-pct_chrs.txt"
    shell:
        r"""
        mkdir -p pri_haps-male/first_hap1/pri pri_haps-male/first_hap1/hap1
        mkdir -p pri_haps-male/first_hap2/pri pri_haps-male/first_hap2/hap2
        cd pri_haps-male
        ln -s ../{input.pri} Mp_RT1Fp.fa
        cd first_hap1
        ln -s ../../{input.hap1} Mh1_RT1Mh2.fa
        cd pri
        {config[repo_location]}/Reciprocal_RagTag.sh ../Mh1_RT1Mh2.fa ../../Mp_RT1Fp.fa {input.reads}
        cd ../hap1
        {config[repo_location]}/Reciprocal_RagTag.sh ../../Mp_RT1Fp.fa ../Mh1_RT1Mh2.fa {input.reads}
        cd ../../first_hap2
        ln -s ../../{input.hap2} Mh2_RT1Mh1.fa
        cd pri
        {config[repo_location]}/Reciprocal_RagTag.sh ../Mh2_RT1Mh1.fa ../../Mp_RT1Fp.fa {input.reads}
        cd ../hap2
        {config[repo_location]}/Reciprocal_RagTag.sh ../../Mp_RT1Fp.fa ../Mh2_RT1Mh1.fa {input.reads}
        cd ../../..
        """

# PASS 3: pri_F (from pass 2 with hap2_F) with hap1_F (pass 1), pri_F (from pass 2 with hap1_F) with hap2_F (pass 1),
#         pri_M (from pass 2 with hap2_M) with hap1_M (pass 1), pri_M (from pass 2 with hap1_M) with hap2_M (pass 1)
rule pass3_F_hap1:
    input:
        pri="pri_haps-female/first_hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        hap="haps_female/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        reads=config["reads"]["female"],
        pDone="pri_haps-female/first_hap2/pri/ragtag_round2-pct_chrs.txt",
        hDone="haps_female/hap1/ragtag_round2-pct_chrs.txt"
    output:
        "pri_haps-female/second-hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-female/second-hap1/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-female/second-hap1/pri/ragtag_round2-pct_chrs.txt",
        "pri_haps-female/second-hap1/hap1/ragtag_round2-pct_chrs.txt"
    shell:
        r"""
        mkdir -p pri_haps-female/second_hap1/pri pri_haps-female/second_hap1/hap1
        cd pri_haps-female/second_hap1
        ln -s ../../{input.pri} Fp_RT2MpFh2.fa
        ln -s ../../{input.hap} Fh1_RT1Fh2.fa
        cd pri
        {config[repo_location]}/Reciprocal_RagTag.sh ../Fh1_RT1Fh2.fa ../Fp_RT2MpFh2.fa {input.reads}
        cd ../hap1
        {config[repo_location]}/Reciprocal_RagTag.sh ../Fp_RT2MpFh2.fa ../Fh1_RT1Fh2.fa {input.reads}
        cd ../../..
        """

rule pass3_F_hap2:
    input:
        pri="pri_haps-female/first_hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        hap="haps_female/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        reads=config["reads"]["female"],
        pDone="pri_haps-female/first_hap1/pri/ragtag_round2-pct_chrs.txt",
        hDone="haps_female/hap2/ragtag_round2-pct_chrs.txt"
    output:
        "pri_haps-female/second-hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-female/second-hap2/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-female/second-hap2/pri/ragtag_round2-pct_chrs.txt",
        "pri_haps-female/second-hap2/hap2/ragtag_round2-pct_chrs.txt"
    shell:
        r"""
        mkdir -p pri_haps-female/second_hap2/pri pri_haps-female/second_hap2/hap2
        cd pri_haps-female/second_hap2
        ln -s ../../{input.pri} Fp_RT2MpFh1.fa
        ln -s ../../{input.hap} Fh2_RT1Fh1.fa
        cd pri
        {config[repo_location]}/Reciprocal_RagTag.sh ../Fh2_RT1Fh1.fa ../Fp_RT2MpFh1.fa {input.reads}
        cd ../hap2
        {config[repo_location]}/Reciprocal_RagTag.sh ../Fp_RT2MpFh1.fa ../Fh2_RT1Fh1.fa {input.reads}
        cd ../../..
        """

rule pass3_M_hap1:
    input:
        pri="pri_haps-male/first_hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        hap="haps_male/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        reads=config["reads"]["male"],
        pDone="pri_haps-male/first_hap2/pri/ragtag_round2-pct_chrs.txt",
        hDone="haps_male/hap1/ragtag_round2-pct_chrs.txt"
    output:
        "pri_haps-male/second-hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-male/second-hap1/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-male/second-hap1/pri/ragtag_round2-pct_chrs.txt",
        "pri_haps-male/second-hap1/hap1/ragtag_round2-pct_chrs.txt"
    shell:
        r"""
        mkdir -p pri_haps-male/second_hap1/pri pri_haps-male/second_hap1/hap1
        cd pri_haps-male/second_hap1
        ln -s ../../{input.pri} Mp_RT2FpMh2.fa
        ln -s ../../{input.hap} Mh1_RT1Mh2.fa
        cd pri
        {config[repo_location]}/Reciprocal_RagTag.sh ../Mh1_RT1Mh2.fa ../Mp_RT2FpMh2.fa {input.reads}
        cd ../hap1
        {config[repo_location]}/Reciprocal_RagTag.sh ../Mp_RT2FpMh2.fa ../Mh1_RT1Mh2.fa {input.reads}
        cd ../../..
        """

rule pass3_M_hap2:
    input:
        pri="pri_haps-male/first_hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        hap="haps_male/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        reads=config["reads"]["male"],
        pDone="pri_haps-male/first_hap1/pri/ragtag_round2-pct_chrs.txt",
        hDone="haps_male/hap2/ragtag_round2-pct_chrs.txt"
    output:
        "pri_haps-male/second-hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-male/second-hap2/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        "pri_haps-male/second-hap2/pri/ragtag_round2-pct_chrs.txt",
        "pri_haps-male/second-hap2/hap2/ragtag_round2-pct_chrs.txt"
    shell:
        r"""
        mkdir -p pri_haps-male/second_hap2/pri pri_haps-male/second_hap2/hap2
        cd pri_haps-male/second_hap2
        ln -s ../../{input.pri} Mp_RT2FpMh1.fa
        ln -s ../../{input.hap} Mh2_RT1Mh1.fa
        cd pri
        {config[repo_location]}/Reciprocal_RagTag.sh ../Mh2_RT1Mh1.fa ../Mp_RT2FpMh1.fa {input.reads}
        cd ../hap1
        {config[repo_location]}/Reciprocal_RagTag.sh ../Mp_RT2FpMh1.fa ../Mh2_RT1Mh1.fa {input.reads}
        cd ../../..
        """

# link final results to output folder. Two options for primaries, prefer version with pass 3 on hap_1 but they should be near identical
rule final:
    input:
        FpA="pri_haps-female/second-hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        Fh1="pri_haps-female/second-hap1/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        fpaDone="pri_haps-female/second-hap1/pri/ragtag_round2-pct_chrs.txt",
        fh1Done="pri_haps-female/second-hap1/hap1/ragtag_round2-pct_chrs.txt",
        FpB="pri_haps-female/second-hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        Fh2="pri_haps-female/second-hap2/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        fpbDone="pri_haps-female/second-hap2/pri/ragtag_round2-pct_chrs.txt",
        fh2Done="pri_haps-female/second-hap2/hap2/ragtag_round2-pct_chrs.txt",
        MpA="pri_haps-male/second-hap1/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        Mh1="pri_haps-male/second-hap1/hap1/ragtag_output-round2/ragtag.scaffold.fasta",
        mpaDone="pri_haps-male/second-hap1/pri/ragtag_round2-pct_chrs.txt",
        mh1Done="pri_haps-male/second-hap1/hap1/ragtag_round2-pct_chrs.txt",
        MpB="pri_haps-male/second-hap2/pri/ragtag_output-round2/ragtag.scaffold.fasta",
        Mh2="pri_haps-male/second-hap2/hap2/ragtag_output-round2/ragtag.scaffold.fasta",
        mpbDone="pri_haps-male/second-hap2/pri/ragtag_round2-pct_chrs.txt",
        mh2Done="pri_haps-male/second-hap2/hap2/ragtag_round2-pct_chrs.txt"       
    output:
        "RT_final/Fp_RTfinal.fa",
        "RT_final/Fh1_RTfinal.fa",
        "RT_final/Fh2_RTfinal.fa",
        "RT_final/Mp_RTfinal.fa",
        "RT_final/Mh1_RTfinal.fa",
        "RT_final/Mh2_RTfinal.fa"
    shell:
        r"""
        mkdir RT_final
        cd RT_final
        ln -s ../{input.FpA} Fp_RTfinal.fa
        ln -s ../{input.Fh1} Fh1_RTfinal.fa
        ln -s ../{input.Fh2} Fh2_RTfinal.fa
        ln -s ../{input.MpA} Mp_RTfinal.fa
        ln -s ../{input.Mh1} Mh1_RTfinal.fa
        ln -s ../{input.Mh2} Mh2_RTfinal.fa
        ln -s ../{input.FpB} ALT-Fp_RTfinal.fa
        ln -s ../{input.MpB} ALT-Mp_RTfinal.fa
        cd ..
        """
