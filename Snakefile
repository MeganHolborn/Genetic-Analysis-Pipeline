import os

# DEFINE CONFIG FILE FOR SNAKEMAKE:
configfile: "config.json"

# DEFINE CONTEXT-VARIABLES:
finalExtensions=['frq.strat', 'fst', 'imiss', 'ld', 'missing.hap', 'hwe', 'het', 'ibc']
locations=set(config['locations'].keys())
samples=set(config['samples'].keys())
superPop=set(config['clusters']['SUPER'])
subPop=set(config['clusters']['SUB'])
bExtensions=["bed", "bim", "fam"]
tExtensions=["map", "ped"]

# Define a list to store a merge-list of dataset names:
mergeList = list()

# def getReferenceGenome(sample):
#     if config['samples'][wildcards.sample]['refGenome'] == "GRCh37" or config['samples'][wildcards.sample]['refGenome'] == "Hg19":
#         return "/nlustre/users/fourie/H.sapiens/gatk_resource_bundle/2.8/hg19/ucsc.hg19.fasta"
#     else if config['samples'][wildcards.sample]['refGenome'] == "GRCh38":
#         return "/nlustre/users/graeme/pipeline-2020/binaries/hg38.fa.gz"


# BEGIN DEFINING RULES:
rule all:
    """
    Catch-all rule to trigger auto-run of all processes. This process will be fired automatically in absence of explicit process name given by cli-argument.
    """
    input:
        expand("final/SUB/ALL_{location}_SUB.{extension}", extension=finalExtensions, location=locations),
        expand("final/SUPER/ALL_{location}_SUPER.{extension}", extension=finalExtensions, location=locations)


# rule VALIDATE:
#     """
#     Perform validation of VCF format as well as REF alleles
#     """

#     input:
#         "input/{{sample}}.vcf.gz"
#     output:

#     params:
#         referenceGenome=getReferenceGenome({wildcards.sample})
#     resources:
#         cpus=28,
#         nodes=1,
#         queue="normal",
#         walltime="30:00:00"
#     run:
#         shell("module load gatk-4.0.12.0; gatk ValidateVariants -R {params.referenceGenome} -V {input}")


rule LIFTOVER:
    """
    Lift Variants onto same Reference build. Otherwise we cant merge them or analyse them in context of eachother.
    """
    input:
        expand("input/{{sample}}.vcf.gz")

    output:
        expand(".intermediates/LIFTOVER/{{sample}}.vcf.gz")
        #expand(".intermediates/LIFTOVER/{{sample}}_LIFTED.unlifted")

    params:
        prefix= lambda wildcards: ".intermediates/LIFTOVER/{sample}_LIFTED".format(sample=wildcards.sample),
        exclusionList= lambda wildcards: ".intermediates/LIFTOVER/{sample}_EXCLUDE.dat".format(sample=wildcards.sample),
        chainFile="binaries/hg19ToHg38.over.chain",
        LiftOver="binaries/liftOverPlink.py",
        rmBadLifts="binaries/rmBadLifts.py",
        ref="binaries/" + config['refGenomes']['GRCh38']

    resources:
        cpus=28,
        nodes=1,
        queue="long",
        walltime="900:00:00"

    run:
        shell("echo 'Determining Liftover requirements now...'")
        if config['samples'][wildcards.sample]['refGenome'] != "GRCh38":
            shell("echo 'Liftover required. All datasets have been mapped to {}'".format(config['samples'][wildcards.sample]['refGenome'])),
            shell("module load liftover"),
            if config['samples'][wildcards.sample]['refGenome'] == "GRCh37" or config['samples'][wildcards.sample]['refGenome'] == "Hg19":
                shell("echo 'Lifting from GRCh37 to GRCh38.'"),
                if not os.path.exists('.intermediates/LIFTOVER'):
                    os.makedirs('.intermediates/LIFTOVER')
                shell("module load plink-2; plink2 --vcf input/{wildcards.sample}.vcf.gz --set-all-var-ids @:#\$r-\$a --allow-extra-chr --new-id-max-allele-len 40 truncate --chr 1-22 --out .intermediates/LIFTOVER/{wildcards.sample}_PREP --export vcf-4.2 bgz --output-chr chr26 --keep-allele-order"),
                shell("sleep 60; tabix -p vcf .intermediates/LIFTOVER/{wildcards.sample}_PREP.vcf.gz"),
                # ToDo: Add in GATK SelectVariant Process and filter out symbolic using `--select-type-to-exclude`:
                shell("module load gatk-4.0.12.0; gatk SelectVariants  -V .intermediates/LIFTOVER/{wildcards.sample}_PREP.vcf.gz --select-type-to-include SNP --select-type-to-include INDEL --select-type-to-exclude MIXED --select-type-to-exclude MNP --select-type-to-exclude SYMBOLIC --exclude-filtered -O .intermediates/LIFTOVER/{wildcards.sample}_CLEANED.vcf.gz"),
                shell("module load picard-2.17.11; java -Xmx128G -jar $PICARD LiftoverVcf I=.intermediates/LIFTOVER/{wildcards.sample}_CLEANED.vcf.gz O=.intermediates/LIFTOVER/{wildcards.sample}.vcf.gz C={params.chainFile} REJECT=.intermediates/LIFTOVER/{wildcards.sample}_REJECTED.vcf.gz R={params.ref}"),
        # TODO: Add conditionals for other human reference genome builds
        else:
            print("No liftover required. Dataset {} is already mapped to GRCh38.".format(wildcards.sample)),
            shell("touch .intermediates/LIFTOVER/{wildcards.sample}_EXCLUDE.dat"),
            shell("module load plink-1.9; plink --map input/{wildcards.sample}.map --ped input/{wildcards.sample}.ped --allow-extra-chr --chr 1-22 --recode vcf --keep-allele-order --exclude {params.exclusionList} --out .intermediates/LIFTOVER/{wildcards.sample}"),
        # shell("bgzip .intermediates/LIFTOVER/{wildcards.sample}.vcf"),
        shell("sleep 1m; tabix -f -p vcf .intermediates/LIFTOVER/{wildcards.sample}.vcf.gz"),
        # shell("echo '.intermediates/LIFTOVER/{wildcards.sample}.vcf.gz' >> .intermediates/LIFTOVER/merge.list")
        mergeList.append(wildcards.sample)



rule ALL_COLLATE:
    """
    Collate Datasets together into 1 psudo-dataset for downstream analysis.
    """
    input:
        expand(".intermediates/LIFTOVER/{sample}.vcf.gz", sample=samples)

    output:
        ".intermediates/COLLATE/ALL.vcf.gz",
        ".intermediates/COLLATE/ALL.vcf.gz.tbi",
        ".intermediates/COLLATE/merge.list"

    params:
        prefix = "ALL_PRE_COLLATE",
        prefix2 = "ALL_MERGE",
        prefix3 = "ALL_SUBFILTERED"

    resources:
        cpus=28,
        nodes=1,
        queue="long",
        walltime="30:00:00"

    run:
        if not os.path.exists('.intermediates/COLLATE'):
            os.makedirs('.intermediates/COLLATE')
        for i in mergeList:
            shell("module load picard-2.17.11; java-jar $PICARD FixVcfHeader I=.intermediates/LIFTOVER/{i}.vcf.gz O=.intermediates/COLLATE/{i}_FIXED.vcf.gz"),
            shell("echo '.intermediates/COLLATE/{i}_FIXED.vcf.gz' >> .intermediates/COLLATE/merge.list"),
        shell("module load bcftools-1.7; bcftools merge -l .intermediates/COLLATE/merge.list -O z -o .intermediates/COLLATE/ALL.vcf.gz"),
        shell("module load samtools-1.7; tabix .intermediates/COLLATE/ALL_INCLUDING_CHR.vcf.gz"),
        shell("module load plink-2; plink2 --vcf .intermediates/COLLATE/ALL_INCLUDING_CHR.vcf.gz --output-chr chr26 --chr 1-22 --expoprt vcf-4.2 bgz --out .intermediates/COLLATE/ALL.vcf.gz")

rule ALL_ANNOTATE:
    """
    Annotate rsID's in psudo-dataset to facilitate down-stream analysis.
    """
    input:
        ".intermediates/COLLATE/ALL.vcf.gz"

    output:
        ".intermediates/ANNOTATE/ALL.vcf.gz",
    
    params:
        refGenome='/apps/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa.gz',
        dbSNP='/nlustre/data/gatk_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz'
    
    resources:
        cpus=28,
        nodes=1,
        queue="normal",
        walltime="30:00:00"

    run:
        # shell("module load plink-1.9; plink --bfile .intermediates/COLLATE/ALL --keep-allele-order --recode --out .intermediates/COLLATE/ALL"),
        shell("module load gatk-4.0.12.0; gatk VariantAnnotator -V {input} -R {params.refGenome} -D /nlustre/data/gatk_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz -O .intermediates/ANNOTATE/ALL_PRE_NAME_CHECKED.vcf"),
        shell("sed -r -e 's/^chr([0-9]{{1,2}})\\t([0-9]+)\\t[0-9]{{1,2}}:[0-9]+[A-Z]{{1}}-[A-Z]{{1}};(rs[0-9]+)/chr\\1\\t\\2\\t\\3/g' .intermediates/ANNOTATE/ALL_PRE_NAME_CHECKED.vcf > .intermediates/ANNOTATE/ALL_ANNOTATED.vcf"),
        shell("module load plink-2; plink2 --vcf .intermediates/ANNOTATE/ALL_ANNOTATED.vcf --export vcf-4.2 bgz --out ALL.vcf"),
        # shell("module load plink-1.9; plink --vcf .intermediates/ANNOTATE/ALL_ANOTATED.vcf --keep-allele-order --make-bed --out .intermediates/ANNOTATE/ALL_ANNOTATED")

# rule Admixture:
#     """
#     Perform Admixture analysis on the large psudo-dataset (Requires 100 000 minimum variants to distinguish sub-populations and 10 000 to distinguish super-populations.)
#     """
#     input:
#         ".intermediates/ANNOTATE/ALL_ANNOTATED.vcf"

#     output:
#         ".intermediates/Admixture/ALL.5.Q",
#         ".intermediates/Admixture/ALL.5.P"

#     params:
#         out_name = ".intermediates/Admixture/ALL",
#         smartPCA = "binaries/EIG-7.2.1/bin/"
    
#     resources:
#         cpus=28,
#         nodes=1,
#         queue="long",
#         walltime="30:00:00"

#     shell:
#         """
#         module load plink-1.9
#         module load admixture-1.3.0
#         plink --vcf {input} --keep-allele-order --thin-count 200000 --set-missing-var-ids @_# --make-bed --out {Params.out_name}
#         admixture {Params.out_name}.bed 5
#         """


rule TRIM_AND_NAME:
    """
    Trim the whole-genome psudo-datasets down to several regions of interest for Variant analysis and Variant effect prediction.
    """
    input:
        ".intermediates/ANNOTATE/ALL.vcf.gz"

    output:
        ".intermediates/TRIM/ALL_{location}_TRIMMED.vcf.gz"

    params:
        fromBP = lambda wildcards: config["locations"][wildcards.location]["GRCh37"]["from"],
        toBP = lambda wildcards: config["locations"][wildcards.location]["GRCh37"]["to"],
        chr = lambda wildcards: config["locations"][wildcards.location]["GRCh37"]["chromosome"]
    
    resources:
        cpus=10,
        nodes=1,
        queue="normal",
        walltime="30:00:00"

    run:
        shell("module load plink-2; plink2 --vcf {input} --from-bp {params.fromBP} --to-bp {params.toBP} --export vcf-4.2 bgz --out .intermediates/TRIM/ALL_{wildcards.location}_TRIMMED"),


rule ALL_FILTER:
    """
    Filter out individuals missing 100% of their variant information (Safety Check).
    """
    input:
        ".intermediates/TRIM/ALL_{location}_TRIMMED.vcf.gz"

    output:
        ".intermediates/FILTER/ALL_{location}_FILTERED.vcf.gz"
    
    resources:
        cpus=10,
        nodes=1,
        queue="short",
        walltime="03:00:00"

    run:
        shell("module load plink-2; plink2 --vcf {input} --mind 1 --output-chr chr26 --export vcf-4.2 bgz--out .intermediates/FILTER/ALL_{wildcards.location}_FILTERED")
        

rule ALL_ANALYZE_SUPER:
    """
    Perform Frequency analysis on super populations.
    """
    input:
        vcf=".intermediates/FILTER/ALL_{location}_FILTERED.vcf.gz",
        popClusters="input/superPopCluster"
    
    output:
        expand("final/SUPER/ALL_{{location}}_SUPER.{extension}", extension=finalExtensions),
        expand("final/SUPER/{i}/ALL_{{location}}_SUPER_{i}_HV.{extensions}", i=superPop, extensions=["log", "ld"])

    params:
        prefix = 'ALL_{location}_SUPER'
    
    resources:
        cpus=15,
        nodes=1,
        queue="normal",
        walltime="30:00:00"

    run:
        shell("module load plink-2; plink2 --vcf {input.vcf} --freq --out final/SUPER/{params.prefix}"),
        shell("module load plink-2; plink2 --vcf {input.vcf} --within {input.popClusters} --freq --fst --missing --indep-pairwise 50 5 .05 --hardy midp --het --out final/SUPER/{params.prefix}"),
        # shell("module load plink-2; plink2 --vcf {input.vcf} --recode HV --out final/SUPER/ALL_{wildcards.location}_SUPER_HV"),
        # shell("module load plink-2; plink2 --vcf {input.vcf} --indep-pairwise 50 10 0.1 --double-id --out final/SUPER/ALL_SUPER_{wildcards.location}"),
        # shell("module load plink-2; plink2 --vcf {input.vcf} --double-id --mind --extract final/SUPER/ALL_SUPER_{wildcards.location}.prune.in --pca var-wts scols=sid --out final/SUPER/ALL_SUPER_{wildcards.location}")
        # for i in superPop:
        #     shell(f"module load plink-1.9; plink --vcf {input.vcf} --double-id --snps-only --keep-allele-order --within {input.popClusters} --keep-cluster-names {i} --r2 inter-chr dprime --recode HV --out final/SUPER/{i}/{params.prefix}_{i}_HV");
        # Admixture:
        # shell("module load plink-1.9; plink --mind --geno --indep-pairwise 50 10 0.1 --vcf {input.vcf} --double-id --keep-allele-order --make-bed --out final/SUPER/ALL_{wildcards.location}_SUPER"),
        # shell("module load admixture-1.3.0; admixture --cv final/SUPER/ALL_{wildcards.location}_SUPER.bed 5"),
        # shell("mv ALL_{wildcards.location}_SUPER.5.* final/SUPER/")


rule ALL_ANALYZE_SUB:
    """
    Perform frequency analysis on sub-populations.
    """
    input:
        vcf=".intermediates/FILTER/ALL_{location}_FILTERED.vcf.gz",
        popClusters="input/subPopCluster"
        
    output:
        expand("final/SUB/ALL_{{location}}_SUB.{extension}", extension=finalExtensions),
        expand("final/SUB/{i}/ALL_{{location}}_SUB_{i}_HV.{extensions}", i=subPop, extensions=["log", "ld"])

    params:
        prefix = 'ALL_{location}_SUB'
    
    resources:
        cpus=15,
        nodes=1,
        queue="normal",
        walltime="30:00:00"
  
    run:
        shell("module load plink-2; plink2 --vcf {input.vcf} --freq --out final/SUB/{params.prefix}"),
        shell("module load plink-2; plink2 --vcf {input.vcf} --within {input.popClusters} --freq --fst --missing --indep-pairwise 50 5 .5 --hardy midp --het --out final/SUB/{params.prefix}"),
        # shell("module load plink-2; plink2 --vcf {input.vcf} --double-id --recode HV --out final/SUB/ALL_{wildcards.location}_SUB_HV")
        # for i in subPop:
        #     shell(f"module load plink-1.9; plink --vcf {input.vcf} --double-id --snps-only --keep-allele-order --within {input.popClusters} --keep-cluster-names {i} --r2 inter-chr dprime --recode HV --out final/SUB/{i}/{params.prefix}_{i}_HV")

# Add in VEP API calls
